/// @cond TYPO

#pragma once

#include <type_traits>
#include <cctype> // isspace

#if __has_include(<string_view>) && defined __cpp_lib_string_view &&\
    __cplusplus >= __cpp_lib_string_view
#   include <string_view>
#else
#   include <sstream>
#endif

template <bool Cond>
using require = typename std::enable_if<Cond>::type;

#if defined __GNUC__ || defined __clang__
#   define PRETTY_FUNC __PRETTY_FUNCTION__
#elif defined _MSC_VER
#   define PRETTY_FUNC __FUNCSIG__
#endif

#ifndef __has_builtin
#   define __has_builtin(x) 0
#endif

#ifndef __has_feature
#   define __has_feature(x) 0
#endif

#ifndef __has_extension
#   define __has_extension(x) 0
#endif

#ifndef __has_cpp_attribute
#   define __has_cpp_attribute(x) 0
#endif

#if __cplusplus/100 >= 2014 || \
   (__has_extension(cxx_generic_lambdas) && __has_extension(cxx_relaxed_constexpr))
#   define TYPO_CPP14_CONSTEXPR_FUNC constexpr
#   define TYPO_ENABLE_CPP14_CONSTEXPR
#else
#   define TYPO_CPP14_CONSTEXPR_FUNC
#endif

#if (defined __cpp_inline_variables && __cplusplus >= __cpp_inline_variables) ||\
    __has_extension(cxx_inline_variables)
#   define TYPO_CPP17_INLINE_VARIABLE inline
#else
#   define TYPO_CPP17_INLINE_VARIABLE
#endif

namespace typo
{

template <typename T>
struct Type;

#if __cplusplus/100 >= 2014
template <typename T>
TYPO_CPP17_INLINE_VARIABLE
constexpr auto type = Type<T> {};
#endif

template <typename T>
struct Type
{
  using type = T;

  constexpr auto remove_cv() const noexcept -> Type<typename std::remove_cv<type>::type>
  {
    return {};
  }
  constexpr auto add_cv() const noexcept -> Type<typename std::add_cv<type>::type>
  {
    return {};
  }

  constexpr auto remove_const() const noexcept -> Type<typename std::remove_const<type>::type>
  {
    return {};
  }
  constexpr auto add_const() const noexcept -> Type<typename std::add_const<type>::type>
  {
    return {};
  }

  constexpr auto remove_volatile() const noexcept -> Type<typename std::remove_volatile<type>::type>
  {
    return {};
  }
  constexpr auto add_volatile() const noexcept -> Type<typename std::add_volatile<type>::type>
  {
    return {};
  }

  constexpr auto remove_reference() const noexcept -> Type<typename std::remove_reference<type>::type>
  {
    return {};
  }
  constexpr auto add_lvalue_reference() const noexcept -> Type<typename std::add_lvalue_reference<type>::type>
  {
    return {};
  }
  constexpr auto add_rvalue_reference() const noexcept -> Type<typename std::add_rvalue_reference<type>::type>
  {
    return {};
  }
  // common shortcuts
  constexpr auto remove_ref() const noexcept -> Type<typename std::remove_reference<type>::type>
  {
    return {};
  }
  constexpr auto add_lvalue_ref() const noexcept -> Type<typename std::add_lvalue_reference<type>::type>
  {
    return {};
  }
  constexpr auto add_rvalue_ref() const noexcept -> Type<typename std::add_rvalue_reference<type>::type>
  {
    return {};
  }

  constexpr auto remove_pointer() const noexcept -> Type<typename std::remove_pointer<type>::type>
  {
    return {};
  }
  constexpr auto add_pointer() const noexcept -> Type<typename std::add_pointer<type>::type>
  {
    return {};
  }
  // common shortcuts
  constexpr auto remove_ptr() const noexcept -> Type<typename std::remove_pointer<type>::type>
  {
    return {};
  }
  constexpr auto add_ptr() const noexcept -> Type<typename std::add_pointer<type>::type>
  {
    return {};
  }

  constexpr auto decay() const noexcept -> Type<typename std::decay<type>::type>
  {
    return {};
  }

  constexpr auto raw() const noexcept -> Type<typename std::remove_reference<typename std::remove_cv<type>::type>::type>
  {
    return {};
  }

  // TODO: Add conversion to string
  operator std::string () const;
};


template <class T1, class T2>
constexpr
bool operator== (Type<T1>, Type<T2>) noexcept
{
  return std::is_same<T1, T2>::value;
}

template <class T1, class T2>
constexpr
bool operator!= (Type<T1> t1, Type<T2> t2) noexcept
{
  return !(t1 == t2);
}



// Printing the Type<T>
namespace detail
{

template <size_t N>
inline
constexpr size_t slen(const char (&) [N])
{
  return N;
}

template <typename T>
inline const char *pretty_name () noexcept
{
  return PRETTY_FUNC;
}


#if __has_include(<string_view>) && defined __cpp_lib_string_view &&\
    __cplusplus >= __cpp_lib_string_view
template <typename T>
#if defined __cpp_consteval && __cplusplus >= __cpp_consteval
consteval
#else
inline constexpr
#endif
auto type_name() -> std::string_view
{
  auto name = std::string_view( pretty_name<T>() );
  auto sq_brace = name.find('[');
  auto T_pos = name.find('T', sq_brace);
  auto eq_sign_pos = name.find('=', T_pos);
  name.remove_prefix(eq_sign_pos + 1);
  name.remove_prefix(std::min(name.find_first_not_of(" "), name.size()));
  name.remove_suffix(name.size() - name.find_first_of("];"));
  return name;
}
#else
#if defined TYPO_ENABLE_CPP14_CONSTEXPR
constexpr size_t pslen(const char *const ps)
{
  for (size_t i = 0; ; ++i)
  {
    if (ps[i] == '\0') return i;
  }
}
#else
constexpr size_t pslen(const char *const ps, size_t idx = 0)
{
  return (ps[idx] == '\0') ?
         idx : pslen(ps, idx + 1);
}
#endif

template <size_t N>
constexpr size_t findChar(const char (&str)[N], char c, size_t idx = 0)
{
  return (str[idx] == c || idx >= N) ?
         idx
         :
         findChar(str, c, idx + 1);

}

template <size_t N>
constexpr size_t skipSpace(const char (&str)[N], size_t idx = 0)
{
  return (!std::isspace( (int)str[idx] )) ?
         idx
         :
         skipSpace(str, idx + 1);
}


struct CSlice
{
  constexpr CSlice (const char *const cstr)
    : from {0}
    , to  {pslen(cstr)}
    , str {cstr}
  {}

  constexpr char operator[] (size_t i) const
  {
    return str[i];
  }

  friend
  std::ostream &operator<< (std::ostream &os, CSlice slice)
  {
    for (size_t i = slice.from; i < slice.to; ++i)
    {
      os << slice[i];
    }

    return os;
  }

  operator std::string() const
  {
    std::stringstream type;

    for (size_t i = from; i < to; i++)
    {
      type << str[i];
    }

    return type.str();
  }

  size_t from;
  size_t to;
  const char *const str;
};

template <typename T>
inline
TYPO_CPP14_CONSTEXPR_FUNC
auto type_name() -> CSlice
{
  auto &&prettyName = PRETTY_FUNC;
  auto name = CSlice( PRETTY_FUNC );
#if defined _MSC_VER
  auto from = findChar(prettyName, '<');
  int angle_brks = 0;
  size_t to = from;

  for (size_t i = from; i < slen(PRETTY_FUNC) - 1; i++)
  {
    if (prettyName[i] == '<') angle_brks += 1;
    else if (prettyName[i] == '>') angle_brks -= 1;

    if (angle_brks == 0)
    {
      to = i;
      break;
    }
  }

  name.from = from + 1;
  name.to = to;
#else
  auto sq_brace = findChar(prettyName, '[');
  auto T_pos = findChar(prettyName, 'T', sq_brace);
  auto eq_sign_pos = findChar(prettyName, '=', T_pos);
  auto type_pos = skipSpace(prettyName, eq_sign_pos + 1);
  name.from = type_pos;// + slen("T = ");
  name.to   = slen(PRETTY_FUNC) - 2;
#endif
  return name;
}
#endif

}//namespace detail


template <typename T>
inline
auto operator<< (std::ostream &os, Type <T>) -> std::ostream &
{
  os << detail::type_name<T>();
  return os;
}

template <typename T>
Type<T>::operator std::string () const
{
  return std::string{ detail::type_name<T>() };
}

#if defined PRETTY_FUNC
#   undef PRETTY_FUNC
#endif
}//namespace typo

/// @endcond
