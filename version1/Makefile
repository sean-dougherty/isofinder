# List of modules in dependency order.
modules=leeadn random special cortes findenes chi_c diver

# The .o files created in bld/ from the modules.
objs=$(addprefix bld/, $(addsuffix .o, $(modules)))

# Create the executable.
isofinder: src/IsoFinderC.f90 $(objs)
	gfortran -Jbld $^ -g -o $@

# Compile a source file.
bld/%.o: src/%.f90 Makefile
	@mkdir -p bld
	gfortran -g -c $< -Jbld -o $@


clean:
	rm -rf bld
	rm -f isofinder
