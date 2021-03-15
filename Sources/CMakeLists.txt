target_sources(xv3rfun
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diagnostic_dot.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/v3rfun_context.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/v3rfun_input.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/v3rfun.f>
)
