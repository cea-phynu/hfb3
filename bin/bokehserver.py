#!/usr/bin/env python3

import os
import numpy as np
import time
import re
import traceback
import select

from bokeh.layouts import row, layout
from bokeh.models import ColumnDataSource, HoverTool, PreText, Div, ColorBar, LinearColorMapper
from bokeh.plotting import figure
from bokeh.colors import RGB
from bokeh.server.server import Server
from bokeh.io import output_file, save
from functools import partial
import numpy as np
import matplotlib.colors

# Path to be created
fifoPath = "/tmp/bokeh.fifo"

from multiprocessing import Process, Manager
import logging

logging.basicConfig(level = logging.INFO, format = '%(asctime)s.%(msecs)03d [%(levelname)s] (%(process)d) %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

maxdata = 100

#===============================================================================
#===============================================================================
#===============================================================================

def bkapp(doc, dataDict, cmdList):

#===============================================================================

  def initPage():
    logging.debug("init page")
    doc.sources = {}
    doc.nbsent = {}
    doc.plots = {}
    doc.theme = 'dark_minimal'
    doc.p0 = True
    del doc.p0
    mainLayout = doc.get_model_by_name('mainLayout')
    t = Div(text = "<h2>Write data to %s</h2>" % (fifoPath))
    mainLayout.children = [row([t], sizing_mode = 'stretch_both')]
    doc.noPlot = True

    cdict = {'red':   [(0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 1.0, 1.0)],
             'green': [(0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 0.0, 0.0)],
             'blue':  [(0.0, 1.0, 1.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 0.0, 0.0)]}

    mpl = matplotlib.colors.LinearSegmentedColormap("tutu", segmentdata = cdict, N = 256)

    palette = mpl(np.linspace(0.0, 1.0, 256))

    doc.bpal = []
    for rgba in palette:
      doc.bpal.append(RGB(r= rgba[0] * 255.0, g= rgba[1] * 255.0, b= rgba[2] * 255.0))

#===============================================================================

  def update():
    cmd = None
    try:
      if doc.lastCmd != len(doc.cmdList):
        cmd = doc.cmdList[doc.lastCmd]
        doc.lastCmd += 1
        logging.debug("============= received command: %s" % (str(cmd)))

        if cmd[0] == 'save': #==================================================
          logging.debug("command: save")
          output_file(cmd[1], mode = 'inline')
          save(doc)
          logging.info("saved in " + cmd[1])

        elif cmd[0] == 'clear': #===============================================
          if cmd[1] == '':
            logging.debug("@clear message")
            logging.info("clear plots")
            initPage()
          else:
            name = cmd[1]
            logging.info("del plot: %s" % (name))

            for r in doc.get_model_by_name('mainLayout').children:
              for p in r.children:
                if p.name == name:
                  logging.info(f"delete {name}")
                  logging.info(f"plots: {r.children}")
                  r.children.remove(p)

                  del doc.sources[name] 
                  logging.info(f"deleted {name}")

        elif cmd[0] == 'test': #================================================
          logging.info("test test test")

        elif cmd[0] == 'updatemap': #===========================================
          name = cmd[1]
          if name in doc.sources:
            logging.debug("command: updatemap")
            newval = dict()
            newval['value'] = [doc.dataDict[name]['z']]

            if 'a' in doc.dataDict[name].keys():
              newval['a'] = [doc.dataDict[name]['a']]

            if 'b' in doc.dataDict[name].keys():
              newval['b'] = [doc.dataDict[name]['b']]

            v = doc.dataDict[name]

            if v['type'] == 'map' or v['type'] == 'mat':

              vmin = np.amin(v['z'])
              vmax = np.amax(v['z'])

              vsym = True
              if abs(vmin) < 1e-14 or abs(vmax) < 1e-14 or vmin * vmax >= 0.0:
                vsym = False

              if vsym:
                if vmax > abs(vmin):
                  vmin = -vmax
                else:
                  vmax = -vmin

              if vsym:
                doc.plots[name].glyph.color_mapper.palette = doc.bpal
                doc.plots[name].glyph.color_mapper.high = vmax
                doc.plots[name].glyph.color_mapper.low = vmin
              else:
                #doc.plots[name].glyph.color_mapper.palette = 'Inferno256'
                doc.plots[name].glyph.color_mapper.palette = 'Turbo256'
                doc.plots[name].glyph.color_mapper.high = vmax
                doc.plots[name].glyph.color_mapper.low = vmin

            doc.sources[name].data = newval

            logging.debug("map source for %s has been updated" % (name))
          else:
            logging.debug("command: updatemap (create)")
            doc.lastCmd -= 1 # retry later

      #=========================================================================

      for k, v in doc.dataDict.items():
        if k not in doc.sources.keys():
          logging.info("new plot: %s" % (k))
          doc.nbsent[k] = 0
          # add plot to page
    
          p = None

          if doc.noPlot:
            doc.get_model_by_name('mainLayout').children = []
            logging.debug("removing message")
            doc.noPlot = False
   
          if v['type'] == 'text': #=============================================
            p = PreText(text = "", width= 200)
            p.name = k
            doc.sources[k] = {'object': p}
            logging.info("new text object %s %d" % (k, v['slot']))

          elif v['type'] == 'map' or v['type'] == 'mat': #========================

            #p = figure(tools = "pan, reset, save, wheel_zoom, box_zoom", active_drag = "pan", active_inspect = None, active_scroll = "wheel_zoom", match_aspect = True)
            p = figure(tools = "pan, reset, save, wheel_zoom, box_zoom", active_drag = "pan", active_scroll = "wheel_zoom", match_aspect = True)
            p.name = k
            
            p.xaxis.axis_label = v['xlabel']
            p.yaxis.axis_label = v['ylabel']
            p.title.text = v['title']

            p.sizing_mode = "stretch_both"
            #p.x_range.range_padding = 0
            #p.y_range.range_padding = 0
            p.toolbar.logo = None
            doc.sources[k] = ColumnDataSource(data = {'value': [v['z']], 'a': [v['a']], 'b': [v['b']]})

            if v['type'] == 'mat':
              p.y_range.start = np.amax(v['a'])
              p.y_range.end   =  0
              p.aspect_scale = 1.0


            cmapper = LinearColorMapper(palette = doc.bpal)

            if v['type'] == 'mat':
              doc.plots[k] = p.image('value', x = v['xmin'], y = v['ymax'], dw = v['xmax'] - v['xmin'], dh = v['ymax'] - v['ymin'], source = doc.sources[k])
            else:
              doc.plots[k] = p.image('value', x = v['xmin'], y = v['ymin'], dw = v['xmax'] - v['xmin'], dh = v['ymax'] - v['ymin'], source = doc.sources[k])

            doc.plots[k].glyph.color_mapper = cmapper

            vmin = np.amin(v['z'])
            vmax = np.amax(v['z'])

            vsym = True
            if abs(vmin) < 1e-14 or abs(vmax) < 1e-14 or vmin * vmax >= 0.0:
              vsym = False

            if vsym:
              if vmax > abs(vmin):
                vmin = -vmax
              else:
                vmax = -vmin

            if vsym:
              doc.plots[k].glyph.color_mapper.palette = doc.bpal
              doc.plots[k].glyph.color_mapper.high = vmax
              doc.plots[k].glyph.color_mapper.low = vmin
            else:
              doc.plots[k].glyph.color_mapper.palette = 'Inferno256'
              doc.plots[k].glyph.color_mapper.high = vmax
              doc.plots[k].glyph.color_mapper.low = vmin

            if v['type'] == 'mat':
              hover_tool = HoverTool(tooltips=[
                                               ('value', '@value'),
                                               ('(a, b)', '(@a, @b)'),
                                              ])
            else:
              hover_tool = HoverTool(tooltips=[
                                               ('value', '@value'),
                                               ('(' + v['xlabel'] + ', ' + v['ylabel'] + ')', '(@b, @a)'),
                                              ])
            p.tools.append(hover_tool)
            p.toolbar.active_inspect = hover_tool

            color_bar = ColorBar(color_mapper=cmapper,
                                 #ticker=BasicTicker(),
                                 label_standoff=10, border_line_color=None, location=(0,0), width=10)

            p.add_layout(color_bar, 'right')

          elif v['type'] == 'temporal': #=======================================
            p = figure(tools = "pan, reset, save, wheel_zoom, box_zoom", active_drag = "pan", active_inspect = None, active_scroll = "wheel_zoom", match_aspect = False)
            p.name = k
            p.x_range.follow = "end"
            p.x_range.range_padding = 0
            p.xaxis.axis_label = v['xlabel']
            p.yaxis.axis_label = v['ylabel']
            p.sizing_mode = "stretch_both"
            p.toolbar.logo = None
            if not hasattr(doc, 'p0'):
              logging.debug("plot is now p0")
              doc.p0 = p
    
            if p is not doc.p0:
              logging.debug("plot is not p0")
              p.x_range = doc.p0.x_range
            doc.sources[k] = ColumnDataSource(data = {'x':[], 'y':[]})
            line   = p.line  (x = 'x', y = 'y', alpha = 1.0, line_width = 1, source = doc.sources[k])
            circle = p.circle(x = 'x', y = 'y', alpha = 1.0, line_width = 1, source = doc.sources[k])
            hover_tool = HoverTool(tooltips=[
                                             ('(x, y)', '(@x, @y)'),
                                            ], renderers=[circle], mode = 'vline')
            p.tools.append(hover_tool)
            p.toolbar.active_inspect = hover_tool
            doc.plots[k] = p

          elif v['type'] == 'xy': #=============================================
            #p = figure(tools = "pan, reset, save, wheel_zoom, box_zoom", active_drag = "pan", active_inspect = None, active_scroll = "wheel_zoom", match_aspect = False)
            p = figure(tools = "pan, reset, save, wheel_zoom, box_zoom", active_drag = "pan", active_inspect = None, active_scroll = "wheel_zoom", match_aspect = False)
            p.name = k
            p.x_range.range_padding = 0
            p.xaxis.axis_label = v['xlabel']
            p.yaxis.axis_label = v['ylabel']
            p.sizing_mode = "stretch_both"
            p.toolbar.logo = None
    
            doc.sources[k] = ColumnDataSource(data = {'x':[], 'y':[]})
            line   = p.line  (x = 'x', y = 'y', alpha = 1.0, line_width = 1, source = doc.sources[k])
            circle = p.circle(x = 'x', y = 'y', alpha = 1.0, line_width = 1, source = doc.sources[k])
            hover_tool = HoverTool(tooltips=[
                                             ('(x, y)', '(@x, @y)'),
                                            ], renderers=[circle], mode = 'vline')
            p.tools.append(hover_tool)
            p.toolbar.active_inspect = hover_tool
            doc.plots[k] = p

          #=====================================================================
          #=====================================================================
          #=====================================================================

          slot = v['slot']

          coords = (0, 0)
          if slot == 0: coords = (0, 0)
          if slot == 1: coords = (0, 1)
          if slot == 2: coords = (1, 0)
          if slot == 3: coords = (1, 1)
          if slot == 4: coords = (2, 0)
          if slot == 5: coords = (2, 1)

          lrows = doc.get_model_by_name('mainLayout').children

          while coords[0] >= len(lrows):
            lrows.append(row(sizing_mode = 'stretch_both'))

          lrow = lrows[coords[0]].children
          while coords[1] >= len(lrow):
            t = PreText(text = "empty", width= 100)
            lrow.append(t)


          lrow[coords[1]] = p

      for k in doc.sources.keys():
        v = {}
        try:
          v = doc.dataDict[k]
        except:
          pass

        if 'type' in v and (v['type'] == 'xy' or v['type'] == 'temporal'):
          vals = v['y'][:]
          if doc.nbsent[k] > len(vals):
            logging.debug("reset data: %s" % (k))
            doc.sources[k].data = {k: [] for k in doc.sources[k].data}
            doc.nbsent[k] = 0

          if doc.nbsent[k] < len(vals):
            nbsent = doc.nbsent[k]
      
            newx = []
            newy = []
      
            logging.debug("sending %d data '%s'" % (len(vals) - nbsent, k))
            for d in vals[nbsent:]:
              newx.append(d[0])
              newy.append(d[1])
      
            new_data = {'x': newx, 'y': newy}
            logging.debug("stream data: %s (%d data)" % (k, len(new_data['x'])))
            if v['type'] == 'temporal':
              doc.sources[k].stream(new_data, maxdata)
            else:
              doc.sources[k].stream(new_data)

            doc.nbsent[k] = len(vals)
        elif 'type' in v and v['type'] == 'text':
          if v['text'] != doc.sources[k]['object'].text:
            doc.sources[k]['object'].text =  v['text']
            logging.info("update text : %s %s" % (k, v['text']))


    except:
      logging.error("in update(): cmd %s", str(cmd))
      logging.error("####################")
      logging.error("#### TRACEBACK #####")
      logging.error("####################")
      logging.error(traceback.print_exc())
      logging.error("####################")
      logging.error("##### CONTINUE #####")
      logging.error("####################")

#===============================================================================

  mainLayout = layout(name = 'mainLayout', sizing_mode = 'stretch_both')
  mainLayout.sizing_mode = 'stretch_both'
  # mainLayout.margin = (-9, -8, -8, -9)
  doc.add_root(mainLayout)

  initPage()

  doc.dataDict = dataDict
  doc.cmdList = cmdList
  doc.lastCmd = 0

  doc.add_periodic_callback(update, 1.0)

#===============================================================================
#===============================================================================
#===============================================================================

def bokeh_process(dataDict, cmdList):
  np.random.seed(1)
  server = Server({'/': partial(bkapp, dataDict = dataDict, cmdList = cmdList)}, num_procs = 1)

  logging.info('Opening Bokeh application on http://localhost:5006/')

  server.start()

  # open local browser
  #server.io_loop.add_callback(server.show, "/", {"browser": 'firefox'})

  server.io_loop.start()

#===============================================================================
#===============================================================================
#===============================================================================

def clearData(d, name):
  dname = "temp_" + name
  if dname in d.keys():
    oldv = d[dname]
    toto = d[dname]['y']
    toto.clear()
    logging.debug("cleared values: " + str(name))
    d[dname]['y'] = toto
    d[dname] = {'slot': oldv['slot'], 'type': oldv['type'], 'y': toto, 'xlabel':oldv['xlabel'], 'ylabel':oldv['ylabel']}
    logging.debug(f"***** d['{dname}'] deleted")
    del(d[dname])


#===============================================================================
#===============================================================================
#===============================================================================

def val(d, xlabel, name, x, y, stype, slot):
  dname = "temp_" + name
  if dname not in d.keys():
    d[dname] = {'slot': slot, 'type': stype, 'y': [], 'xlabel':xlabel, 'ylabel':name}
  toto = d[dname]['y']
  toto.append((x, y))
  logging.debug("values: " + str(toto))
  d[dname] = {'slot': slot, 'type': stype, 'y': toto, 'xlabel':xlabel, 'ylabel':name}

#===============================================================================
#===============================================================================
#===============================================================================

def line_split(line):
  return re.findall(r'[^"\s]\S*|".+?"', line)

#===============================================================================
#===============================================================================
#===============================================================================

if __name__ == '__main__':
  with Manager() as manager:
    currentMap = None
    dataDict = manager.dict()
    cmdList = manager.list()

    bokeh_process = Process(target = bokeh_process, args = (dataDict, cmdList)).start()

    try:
      os.mkfifo(fifoPath)
    except:
      logging.info("Failed to create FIFO")

    t0 = time.time()
    currentName = None
    currentSlot = -1
    oldName = None

    #===========================================================================

    with open(fifoPath) as fifo:
      while True:
        select.select([fifo],[],[fifo])
        for mesg in fifo:
          try:
            words = line_split(mesg[:-1])
            logging.debug("words: %s" % (str(words)))
            cmd = words[0]
            if cmd[0] == '@':
              cmd = cmd[1:]
              if cmd == "clear": #===============================================
                name = ''
                if (len(words) > 1):
                  name = words[1].strip('"')

                logging.debug("===== name: " + name)
                if name == '':
                  dataDict.clear()
                  cmdList.append(("clear", name))
                  t0 = time.time()
                  currentName = None
                  currentSlot = -1
                  oldName = None
                else:
                  clearData(dataDict, name)
                  # dname = "temp_" + name
                  # cmdList.append(("clear", dname))

              elif cmd == "slot": #==============================================
                if (len(words) > 1):
                  currentSlot = int(words[1]) - 1
                  oldName = None

              elif cmd == "text": #==============================================
                text = ''
                if (len(words) > 2):
                  currentName = words[1].strip('"')

                  if oldName != currentName:
                    oldName = currentName
                    currentSlot += 1
                    logging.debug("slot: %d name: %s" % (currentSlot, currentName))
  
                  text = words[2].strip('"')
                  if currentName in dataDict:
                    if 'text' in dataDict[currentName]:
                      currentText = dataDict[currentName]['text'][:]
                      # text = currentText + text

                  dataDict[currentName] = {'slot': currentSlot, 'type': 'text', 'text': text}

              elif cmd == "save": #==============================================
                name = 'result.html'
                if (len(words) > 1):
                  name = words[1].strip('"')
                cmdList.append(("save", name))

              elif cmd == "xy_data": #===========================================
                logging.debug("words: %r" % (words))
                xlabel = words[1].strip('"')
                ylabel = words[2].strip('"')
                x = float(words[3])
                y = float(words[4])

                currentName = ylabel

                if oldName != currentName:
                  oldName = currentName
                  currentSlot += 1
                  logging.debug("slot: %d name: %s" % (currentSlot, currentName))

                val(dataDict, xlabel, currentName, x, y, "xy", currentSlot)
              elif cmd == "temp_data": #=========================================
                if len(words) == 4:
                  x = float(words[2])
                  y = float(words[3])
                else:  
                  x = time.time() - t0
                  y = float(words[2])
                ylabel = words[1].strip('"')

                currentName = ylabel

                if oldName != currentName:
                  currentSlot += 1
                  oldName = currentName
                  logging.debug("slot: %d name: %s" % (currentSlot, currentName))

                val(dataDict, "time [s]", currentName, x, y, "temporal", currentSlot)
              elif cmd == "mat": #===============================================
                name = words[1].strip('"')
                nbx  = int(words[2])
                nby  = int(words[3])
                iy   = int(words[4])

                currentName = 'mat_' + name

                if oldName != currentName:
                  oldName = currentName
                  currentSlot += 1
                  logging.debug("slot: %d name: %s" % (currentSlot, currentName))

                iw = 5
                if iy == 0:
                  currentMap = np.zeros((nbx, nby), dtype = np.float32)
                  aval = np.zeros((nbx, nby), dtype = np.int16)
                  bval = np.zeros((nbx, nby), dtype = np.int16)
  
                for ix in range(nbx):
                  currentMap[(nbx - 1 - ix, iy)] = float(words[iw])
                  aval[(ix, iy)] = nbx - 1 - ix
                  bval[(ix, iy)] = iy
                  iw += 1
  
                if iy == nby - 1:
                  xmin = -0.5
                  ymin = -0.5
                  xmax = nbx - 1 + 0.5
                  ymax = nby - 1 + 0.5
                  dataDict[currentName] = {'slot': currentSlot, 'type': 'mat', 'z': currentMap, 'a': aval, 'b': bval, 'xmin': xmin, 'ymin': ymin, 'xmax': xmax, 'ymax': ymax, 'title': name, 'xlabel': 'b', 'ylabel': 'a'}
                  cmdList.append(("updatemap", currentName))
              elif cmd == "map": #===============================================
                name = words[1].strip('"')
                xlabel = words[2].strip('"')
                ylabel = words[3].strip('"')
                nbx  = int(words[8])
                nby  = int(words[9])
                iy   = int(words[10])
                iw = 11

                currentName = 'map_' + name

                if oldName != currentName:
                  oldName = currentName
                  currentSlot += 1
                  logging.debug("slot: %d name: %s" % (currentSlot, currentName))

                if currentMap is None or iy == 0:
                  currentMap = np.zeros((nbx, nby), dtype = np.float32)
                  aval = np.zeros((nbx, nby))
                  bval = np.zeros((nbx, nby))
  
                for ix in range(nbx):
                  currentMap[(ix, iy)] = float(words[iw])
                  iw += 1
                  bval[(ix, iy)] = ix
                  aval[(ix, iy)] = iy
  
                if iy == nby - 1:
                  xmin = float(words[4])
                  ymin = float(words[5])
                  xmax = float(words[6])
                  ymax = float(words[7])
                  stepx = 1.0
                  stepy = 1.0
                  if nbx != 1: stepx = (xmax - xmin) / (float(nbx) - 1.0)
                  if nby != 1: stepy = (ymax - ymin) / (float(nby) - 1.0)

                  bval = (bval * stepx) + xmin
                  aval = (aval * stepy) + ymin

                  xmin -= stepx / 2.0
                  ymin -= stepy / 2.0
                  xmax += stepx / 2.0
                  ymax += stepy / 2.0

                  currentMap = np.transpose(currentMap)
                  aval = np.transpose(aval)
                  bval = np.transpose(bval)

                  dataDict[currentName] = {'slot': currentSlot, 'type': 'map', 'z': currentMap, 'a': aval, 'b': bval, 'xmin': xmin, 'ymin': ymin, 'xmax': xmax, 'ymax': ymax, 'title': name, 'xlabel': xlabel, 'ylabel': ylabel}

                  cmdList.append(("updatemap", currentName))

              elif cmd == "test": #==============================================
                  cmdList.append(("test",))

            else:
              currentName = 'noname'

              if oldName != currentName:
                oldName = currentName
                currentSlot += 1
                logging.debug("slot: %d name: %s" % (currentSlot, currentName))

              if len(words) == 2:
                x = float(words[0])
                y = float(words[1])
                val(dataDict, "X axis", currentName, x, y, "xy", currentSlot)
              else:  
                x = time.time() - t0
                y = float(words[0])
                val(dataDict, "time [s]", currentName, x, y, "temporal", currentSlot)
          except:
            logging.error("Bad input: %s", mesg[:-1])
            logging.error("####################")
            logging.error("#### TRACEBACK #####")
            logging.error("####################")
            logging.error(traceback.print_exc())
            logging.error("####################")
            logging.error("##### CONTINUE #####")
            logging.error("####################")
            #sys.exit(0)

  bokeh_process.join()
