make clean
make
./main inputs/Silica_beta_cristobalite.dat output/data/MadelungConstant_silica.dat
./main inputs/TiO2_Rutile.dat output/data/MadelungConstant_rutile.dat
./main inputs/ZnS_Sphalerite.dat output/data/MadelungConstant_ZnS.dat
make -C output/gnuplot

