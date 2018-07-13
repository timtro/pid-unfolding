macro (create_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro (create_test)

set(testname pid-frp)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH "-lsodium")
    if (WITH_PLOT)
    target_link_libraries (${testname} "-lboost_iostreams" "-lboost_system"
      "-lboost_filesystem")
    endif()
    create_test(${testname})

set(testname pid-unfold)
    add_executable (${testname} ${testname}.cpp)
    target_link_libraries (${testname} CATCH ${RANGE-V3_LIBRARIES} "-lsodium")
    if (WITH_PLOT)
    target_link_libraries (${testname} "-lboost_iostreams" "-lboost_system"
      "-lboost_filesystem")
    endif()
    create_test(${testname})
