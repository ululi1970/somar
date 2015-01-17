# An executable will be created for each name in ebase.
ebase := somar

# Where do the Chombo libraries live?
CHOMBO_HOME = ../../lib/Chombo/lib

# Do you want to trap floating point errors?
TRAPFPE = FALSE


# BEGIN: Don't mess with this stuff.
MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config
LibNames :=  AMRElliptic AMRTools AMRTimeDependent BoxTools

base_dir =  .
src_dirs = $(shell find ../src -maxdepth 256 -type d | tr '\n' ' ')
src_dirs+= ../DNS/projects/Beam_Profile
inc_dirs = $(CHOMBO_HOME)/include $(src_dirs) /home/eds/research/somar/DNS/src

include $(CHOMBO_HOME)/mk/Make.example
# END: Don't mess with this stuff.
