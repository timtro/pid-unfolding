macro (create_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro (create_test)

set(testname utest-numerics-pid_algebra)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})

set(testname utest-numerics-Plant)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})

set(testname utest-unfold)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH)
    create_test(${testname})
