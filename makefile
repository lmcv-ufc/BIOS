############################################################################
#                                                            
# Makefile do programa BIOS
#                                                            
############################################################################

# Application constants

DIR = .
CTRL		= $(DIR)/ctrl
MAIN		= $(DIR)/main
MLIB		= $(DIR)/mlib
OBJ             = $(DIR)/obj
PEN		= $(DIR)/pen
PROB		= $(DIR)/prob
SOL		= $(DIR)/sol
SUR		= $(DIR)/sur

# Directory constants

DIRS =	$(CTRL) $(PROB) $(MAIN) $(MLIB) $(PEN) \
        $(SOL) $(SUR)

# Compilation parameters

#CC	= mpiCC
CC  = g++ 
POLICE	= -Wall
DEBUG	= 
DEBUG	= -g
#DEBUG	= -g -pg
OPTMIZ	=
OPTMIZ	= -O2
#DEFINE	= -D_UNIX_ -D_MPI_ -D_OMP_
#DEFINE = -D_UNIX_ -D_OMP_
DEFINE	= -D_UNIX_
INCLUDE	= -I$(CTRL) -I$(PROB) -I$(MAIN) -I$(MLIB) -I$(PEN) \
          -I$(SOL) -I$(SUR)
#CFLAGS	= -fopenmp $(INCLUDE) $(POLICE) $(DEBUG) $(OPTMIZ) $(DEFINE)
CFLAGS	= $(INCLUDE) $(POLICE) $(DEBUG) $(OPTMIZ) $(DEFINE)
SYSLIBS	= -lm  

# Aplication modules

CTRLMOD	=		\
	optalg		\
	stdga           \
	stdpso          \
	stdabc          \
	stdais          \
	lamalg		\
	lamga		\
	algparam	\
	modpso          \
	sel		\
	sao		\
	saorbf		\
	saokrg		\
	modnsgaII	\
	modlamnsgaII	\
	stdde		\
	rs

	#kitayamasao	\
	#grego		\


INDMOD =		\
	optsolution	\
	individual	\
	particle	\
	food		\
	sampsao		\
	group

MAINMOD	=		\
	input		\
	main		\
	utl

MLIBMOD	=		\
	mat		\
	matvec		\
	vec		\
	rk4			

PENMOD	=		\
	penalty

SURMOD	=		\
	krg		\
	rbf		\
	surr 		\
	samp		\
	problike	

PROBMOD	=		\
	benchmark	\
	lam	        \
	lamplt		\
	fgm		\
	fgmplt		\
	material 	\
	metaopt 	\
	metaprob 	\
	problem		\
        probsurr	


# Object modules

OBJS	= $(CTRLMOD:%=$(OBJ)/%.o)	\
	  $(INDMOD:%=$(OBJ)/%.o)	\
	  $(MAINMOD:%=$(OBJ)/%.o)	\
	  $(MLIBMOD:%=$(OBJ)/%.o)	\
	  $(PENMOD:%=$(OBJ)/%.o)	\
	  $(SURMOD:%=$(OBJ)/%.o)	\
	  $(PROBMOD:%=$(OBJ)/%.o)

# Executable (To use the profile, use '-static -pg' )

FIRST: $(DIRS)
	$(MAKE) $(DIR)/bios

$(DIR)/bios	: $(DIRS) $(OBJS)
		$(CC) -o $@ $(OBJS) $(SYSLIBS)
#	$(CC) -fopenmp -o $@ $(OBJS) $(SYSLIBS)
#		$(CC) -o $@ $(OBJS) -static -pg $(SYSLIBS)

# Object compilation

$(OBJ)/%.o	: $(CTRL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MAIN)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MLIB)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(PEN)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(SUR)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(PROB)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(SOL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

# Dependencies generation

depend:	FORCE
	@echo ""
	@echo "Generating CTRL dependencies..."
	-@g++ $(INCLUDE) -MM $(CTRL)/*.cpp > biosdep.tmp
	@echo "Generating MAIN dependencies..."
	-@g++ $(INCLUDE) -MM $(MAIN)/*.cpp >> biosdep.tmp
	@echo "Generating MLIB dependencies..."
	-@g++ $(INCLUDE) -MM $(MLIB)/*.cpp >> biosdep.tmp
	@echo "Generating PEN dependencies..."
	-@g++ $(INCLUDE) -MM $(PEN)/*.cpp >> biosdep.tmp
	@echo "Generating SUR dependencies..."
	-@g++ $(INCLUDE) -MM $(SUR)/*.cpp >> biosdep.tmp
	@echo "Generating PROB dependencies..."
	-@g++ $(INCLUDE) -MM $(PROB)/*.cpp >> biosdep.tmp
	@echo "Generating SOL dependencies..."
	-@g++ $(INCLUDE) -MM $(SOL)/*.cpp >> biosdep.tmp
	@sed -e '1,$$s/^\([^ ]\)/$$(OBJ)\/\1/' < biosdep.tmp > biosdep && \
          rm -f biosdep.tmp
	@echo "Generation completed."

FORCE:

# dependencies


###################################################### End of file #########
