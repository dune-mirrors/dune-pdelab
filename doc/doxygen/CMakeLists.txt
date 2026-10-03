# Create and install doxygen documentation
add_doxygen_target(THEME awesome
  DEPENDS Mainpage.md modules.txt bibliography.bib)
