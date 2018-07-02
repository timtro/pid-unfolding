macro (create_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro (create_test)

add_library (CATCH catch.cpp)

find_package ( Boost REQUIRED )
include_directories ( ${Boost_INCLUDE_DIRS} )

set(testname utest-numerics-pid_algebra)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})

set(testname utest-numerics-Plant)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})

set(testname utest-logic-frp)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH "-lsodium")
    create_test(${testname})

set(testname utest-pid-frp)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH "-lsodium")
    if (WITH_PLOT)
    target_link_libraries (${testname} "-lboost_iostreams" "-lboost_system"
      "-lboost_filesystem")
    endif()
    create_test(${testname})
