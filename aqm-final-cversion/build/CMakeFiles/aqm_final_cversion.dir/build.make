# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gk141/Documents/Projects/PHY765/aqm-final-cversion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build

# Include any dependencies generated for this target.
include CMakeFiles/aqm_final_cversion.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/aqm_final_cversion.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/aqm_final_cversion.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/aqm_final_cversion.dir/flags.make

CMakeFiles/aqm_final_cversion.dir/main.cpp.o: CMakeFiles/aqm_final_cversion.dir/flags.make
CMakeFiles/aqm_final_cversion.dir/main.cpp.o: /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/main.cpp
CMakeFiles/aqm_final_cversion.dir/main.cpp.o: CMakeFiles/aqm_final_cversion.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/aqm_final_cversion.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/aqm_final_cversion.dir/main.cpp.o -MF CMakeFiles/aqm_final_cversion.dir/main.cpp.o.d -o CMakeFiles/aqm_final_cversion.dir/main.cpp.o -c /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/main.cpp

CMakeFiles/aqm_final_cversion.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/aqm_final_cversion.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/main.cpp > CMakeFiles/aqm_final_cversion.dir/main.cpp.i

CMakeFiles/aqm_final_cversion.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/aqm_final_cversion.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/main.cpp -o CMakeFiles/aqm_final_cversion.dir/main.cpp.s

# Object files for target aqm_final_cversion
aqm_final_cversion_OBJECTS = \
"CMakeFiles/aqm_final_cversion.dir/main.cpp.o"

# External object files for target aqm_final_cversion
aqm_final_cversion_EXTERNAL_OBJECTS =

aqm_final_cversion: CMakeFiles/aqm_final_cversion.dir/main.cpp.o
aqm_final_cversion: CMakeFiles/aqm_final_cversion.dir/build.make
aqm_final_cversion: CMakeFiles/aqm_final_cversion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable aqm_final_cversion"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/aqm_final_cversion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/aqm_final_cversion.dir/build: aqm_final_cversion
.PHONY : CMakeFiles/aqm_final_cversion.dir/build

CMakeFiles/aqm_final_cversion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/aqm_final_cversion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/aqm_final_cversion.dir/clean

CMakeFiles/aqm_final_cversion.dir/depend:
	cd /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gk141/Documents/Projects/PHY765/aqm-final-cversion /home/gk141/Documents/Projects/PHY765/aqm-final-cversion /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build /home/gk141/Documents/Projects/PHY765/aqm-final-cversion/build/CMakeFiles/aqm_final_cversion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/aqm_final_cversion.dir/depend
