macro (create_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro (create_test)

add_library (CATCH catch.cpp)

find_package ( Boost REQUIRED )
include_directories ( ${Boost_INCLUDE_DIRS} )

set(testname utest-pid-simple)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})

set(testname utest-pid-simulations)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})
