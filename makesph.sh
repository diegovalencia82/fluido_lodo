rm *.o
rm *.a
rm sph


gfortran -c input.f
gfortran -c time_integration.f
gfortran -c single_step.f
gfortran -c neighboring_search.f
gfortran -c kernel.f
gfortran -c wijdwij.f
gfortran -c density.f
gfortran -c presioni.f
gfortran -c momento.f
gfortran -c sph_presion.f
gfortran -c tasa_deformacion_epsilon.f
gfortran -c viscous_force.f
gfortran -c external_force.f
gfortran -c indexx.f
gfortran -c neighboring_searchx.f

gfortran -c tasa_deformacion_epsilon2.f
gfortran -c viscous_force_mud.f
gfortran -c viscous_bingham.f


ar ru mysph.a input.o time_integration.o single_step.o neighboring_search.o kernel.o wijdwij.o density.o presioni.o momento.o sph_presion.o tasa_deformacion_epsilon.o viscous_force.o external_force.o indexx.o neighboring_searchx.o tasa_deformacion_epsilon2.o  viscous_force_mud.o viscous_bingham.o

make sph
