###########################################################
# Rule
###########################################################
# target: prerequisites
#	command
###########################################################
# Phony targets are not real files, they are simply names for actions to be performed
###########################################################
#    @ 	In front of a command prevents the printing of the command itself
#    $@ Represents the target of the rule.
#    $^ Represents all prerequisites of the rule, without duplicates.
#    $< Represents the first prerequisite of the rule.
#    $? Represents all prerequisites newer than the target.
#    $  Represents the stem of the target or prerequisite when using pattern rules.
#    := Evaluation when the variable in declared
#    =  Evaluation when the variable in called
###########################################################
#                       To edit:
###########################################################

#compiler
CXX := g++
#flags
CXXFLAGS := -std=c++17 -Wall -Wextra -MMD -fopenmp -g
#Executable
EXEC := hydrogen.x
#Source files
SOURCE_DIRS := src lib/src
#Header file
INCLUDE_DIRS := lib/include
###########################################################
#          	     Used variables:
###########################################################
#Objects and dependencies directory (invisible directories to be created inside the src directories)
BUILD_DIR := .build
BUILD_DIRS := $(foreach src_dir, $(SOURCE_DIRS), $(src_dir)/$(BUILD_DIR))
SOURCE_FILES := $(foreach src_dir, $(SOURCE_DIRS), $(wildcard $(src_dir)/*.cpp) )                            
INCLUDES := $(foreach incl, $(INCLUDE_DIRS), -I$(incl) )                     
OBJECT_FILES := $(foreach src_file, $(SOURCE_FILES), $(dir $(src_file))$(BUILD_DIR)/$(notdir $(src_file)).o)
DEPENDENCY_FILES := $(OBJECT_FILES:.o=.d)

###########################################################
#  		      Build Rules:
###########################################################
#all target:
.PHONY: all
all: $(BUILD_DIRS) $(EXEC)
	@echo "********************************"
	@echo "********* All compiled *********"
	@echo "********************************"
#----------------------------------------------------------
# Build target:
$(EXEC): $(OBJECT_FILES)
	$(CXX) $^ -o $@ $(CXXFLAGS) $(INCLUDES)
#----------------------------------------------------------
# include all Dependencies:
-include $(DEPENDENCY_FILES)

#this macro takes only one argument with a .cpp file
#Here we must use $$< and $$@ to escape variables $< and $@
define BUILD_OBJECT_RULES
$(dir $(1))$(BUILD_DIR)/$(notdir $(1)).o : $(1) | $(BUILD_DIRS)
	$(CXX) -c $$< -o $$@ $(CXXFLAGS) $(INCLUDES)
endef
#.o rules:
$(foreach file,$(SOURCE_FILES),$(eval $(call BUILD_OBJECT_RULES,$(file))))
#----------------------------------------------------------
$(BUILD_DIRS):
	mkdir -p $(BUILD_DIRS)
#----------------------------------------------------------
.PHONY: clean
clean:
	rm -rf $(BUILD_DIRS) $(EXEC)

.PHONY: run
run: all
	@echo "********************************"
	@echo "*********** Program ************"
	@echo "********************************"
	@./$(EXEC)