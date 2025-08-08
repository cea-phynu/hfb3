#!/usr/bin/env python3

# ==============================================================================
# ==============================================================================
# ==============================================================================

import csv

# ==============================================================================
# ==============================================================================
# ==============================================================================


def gen_agr(resultFile, template, files, columnLabel):

    with open(resultFile, 'w') as result:
        headerStr = open(template, 'r').read()

        result.write(headerStr)

        for id, fileName in enumerate(files):
            result.write(f"@target G0.S{id}\n")
            result.write("@type xy\n")
            with open(fileName, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    result.write(f"{row['beta2']} {row[columnLabel]}\n")
            result.write("&\n")
        result.close()

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == '__main__':

    filesBasis = ('240Pu_1x14_Q1.00.csv',
                  '240Pu_1x16_Q1.00.csv',
                  '240Pu_1x18_Q1.00.csv',
                  '240Pu_2x11_Q1.00.csv',
                  '240Pu_2x12_Q1.00.csv',
                  '240Pu_2x13_Q1.00.csv',
                  )

    # figure 3
    gen_agr('240Pu_relative.agr'     , '240Pu_relative.agr.header'     , filesBasis, 'Ehfb normalized')

    # figure 4
    gen_agr('240Pu_relative_zoom.agr', '240Pu_relative_zoom.agr.header', filesBasis, 'Ehfb normalized')

    # figure 5
    gen_agr('240Pu_absolute.agr'     , '240Pu_absolute.agr.header'     , filesBasis, 'Ehfb')

    # figure 6
    gen_agr('240Pu_absolute_zoom.agr', '240Pu_absolute_zoom.agr.header', filesBasis, 'Ehfb')

    filesGQ = ('240Pu_1x16_Q1.00.csv',
               '240Pu_1x16_Q1.10.csv',
               '240Pu_1x16_Q1.20.csv',
               '240Pu_1x16_Q1.30.csv',
               '240Pu_2x11_Q1.00.csv',
               '240Pu_2x11_Q1.30.csv',
               )

    # figure 7
    gen_agr('240Pu_gq.agr', '240Pu_gq.agr.header', filesGQ, 'Ehfb normalized')
