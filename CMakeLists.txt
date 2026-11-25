cmake_minimum_required(VERSION 3.12)
project(Coronarias)

find_package(VTK REQUIRED COMPONENTS
  CommonCore CommonDataModel CommonMath CommonExecutionModel
  IOImage IOLegacy
  FiltersCore FiltersGeneral FiltersModeling FiltersSources
  InteractionStyle RenderingCore RenderingOpenGL2 RenderingFreeType)

add_executable(coronarias main.cpp)
target_link_libraries(coronarias PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS coronarias MODULES ${VTK_LIBRARIES})
