if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(PEDANTIC_COMPILE_FLAGS
      -pedantic-errors
      -Wall
      -Wextra
      -pedantic
      -Wold-style-cast
      -Wundef
      -Wredundant-decls
      -Wwrite-strings
      -Wpointer-arith
      -Wcast-qual
      -Wmissing-include-dirs
      -Wcast-align
      -Wctor-dtor-privacy
      -Wdisabled-optimization
      -Winvalid-pch
      -Woverloaded-virtual
      -Wconversion
      -Wundef
      -Wno-ctor-dtor-privacy
      -Wno-format-nonliteral)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(PEDANTIC_COMPILE_FLAGS
      -Wall
      -Werror
      -Wextra
      -pedantic
      -Wconversion
      -Wundef
      -Wdeprecated
      -Wweak-vtables
      -Wshadow
      -Wno-gnu-zero-variadic-macro-arguments)
endif()

set(OPTIM_COMPILE_FLAGS -O3 -march=native -ffast-math)
