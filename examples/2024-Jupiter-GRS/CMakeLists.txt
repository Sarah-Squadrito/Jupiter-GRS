# ===========================================
# Example Problem of Simulating Jupiter's GRS
# ===========================================

# 1. Compile jupiter crm
setup_problem(jupiter_grs)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
