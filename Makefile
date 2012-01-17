C_LIBS=-lGL ./libglfw3.a
C_CPP_COMMON_COMPILE_FLAGS:= -O3 -g -Wall -Wextra -Wuninitialized -Winit-self -Wstrict-aliasing -Wfloat-equal -Wshadow -Wconversion -Werror -fpack-struct=4

C_COMPILE:=gcc -c
C_COMPILE_FLAGS:=-ansi -pedantic-errors

C_LINK:=g++
C_LINK_FLAGS:=-g -lm -lpthread $(C_LIBS)

DECODE_DEPENDS:=window.c
C_PROJECTS:=decode

all: $(C_PROJECTS)

define upperString
$(shell echo $1 | tr [a-z] [A-Z] )
endef

define C_PROJECT_template
$2_SRCFILES += $1.c
$2_SRCFILES += $($2_DEPENDS)
$2_DEPEND_OBJS:=$($2_DEPENDS:.c=.o)

$2_OBJFILE:=$1.o
$2_OBJFILES:=$$($2_SRCFILES:.c=.o)

C_SRCFILES += $$($2_SRCFILES)
C_OBJFILES += $$($2_OBJFILES)

C_TARGETS += $1

$$($2_OBJFILE): $$($2_DEPEND_OBJS) $1.c
$1: $$($2_OBJFILES) 
endef
     
$(foreach project,$(C_PROJECTS),$(eval $(call C_PROJECT_template,$(project),$(call upperString,$(project)))))

test:
	@echo C_PROJECTS=$(C_PROJECTS)
	@echo C_TARGETS=$(C_TARGETS)
	@echo C_SRCFILES=$(C_SRCFILES)
	@echo C_OBJFILES=$(C_OBJFILES)

%.8: %.go
	@echo Go Compiling $<
	@$(GC) $<

%: %.8
	@echo Go Linking $<
	@$(LD) -o $@ $<

%.o: %.c
	@echo C Compiling $<
	@$(C_COMPILE) -MMD $(C_CPP_COMMON_COMPILE_FLAGS) $(C_COMPILE_FLAGS) -o $*.o $<

%: %.o
	@echo C Linking $<
	@$(C_LINK) -o $@ $^ $(C_LINK_FLAGS)

.PHONY: all clean nuke 
.SUFFIXES:            # Delete the default suffixes

FORCE:

clean: FORCE
	rm -f $(C_OBJFILES)

nuke: clean
	rm -f $(C_TARGETS)
