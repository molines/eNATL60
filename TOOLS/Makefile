# make file for tools

include make.macro

ALL: mask_old_zone.exe

mask_old_zone.exe: mask_old_zone.f90
	$(F90) $(FFLAGS) mask_old_zone.f90 -o  mask_old_zone.exe
