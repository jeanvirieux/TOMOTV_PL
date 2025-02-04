#
#
donothing: # To protect against an idle "make" to "see what happens"
	@echo "Type make all for compiling the different modules or"
	@echo "make model for the module MODEL"
	@echo "make raytrace for the module RAYTRACING"
	@echo "make derivation for the module DERIVE"
	@echo "make precond for the module PRECOND"
	@echo "make inversion for the module INVERSE"
	@echo "....................................."
	@echo "make inversion for the module UTIL"
	@echo "make inversion for the module STATISTICS"
	@echo "====================================="

all : modelisation raytrace derivation precondition inversion utilities stats thirdparties

modelisation : 
	cd ./MODEL; $(MAKE)

raytrace : 
	cd ./RAYTRACING; $(MAKE)

derivation : 
	cd ./DERIVE; $(MAKE)

precondition : 
	cd ./PRECOND; $(MAKE)

inversion : 
	cd ./INVERSE; $(MAKE)

utilities : 
	cd ./UTIL; $(MAKE) 

stats : 
	cd ./STATISTICS; $(MAKE) 

thirdparties : 
	cd ./THIRD_PARTY; $(MAKE) 


clean :
	cd ./MODEL; $(MAKE) clean
	cd ./RAYTRACING; $(MAKE) clean
	cd ./DERIVE; $(MAKE) clean
	cd ./PRECOND; $(MAKE) clean
	cd ./INVERSE; $(MAKE) clean
	cd ./UTIL; $(MAKE) clean
	cd ./STATISTICS; $(MAKE) clean
	cd ./THIRD_PARTY; $(MAKE) clean
	rm -rf ../BIN_LINUX_PL/*

install :
	cd ./MODEL; $(MAKE) clean
	cd ./RAYTRACING; $(MAKE) clean
	cd ./DERIVE; $(MAKE) clean
	cd ./PRECOND; $(MAKE) clean
	cd ./INVERSE; $(MAKE) clean
	cd ./UTIL; $(MAKE) clean
	cd ./STATISTICS; $(MAKE) clean
	cd ./THIRD_PARTY; $(MAKE) clean

