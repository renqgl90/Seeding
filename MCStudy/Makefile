SRC = src
LIB = lib
BIN = bin
INC = $(addprefix -I,$(SRC))

vpath %.C $(SRC)

CC      = g++
CCFLAGS = -c  -Wno-deprecated -g

LD      = g++
LDFLAGS = 

ROOTINCS = $(shell root-config --cflags) # will be used as it is
ROOTLIBS = $(shell root-config --glibs)  # will be used as it is
LIBS     = $(LIB)                        # will get an -L for each entry

EXEC0        = RunSelectorXLayers
COMPONENTS0  = RunSelectorXLayers FakedSeeding #FakedSeeding #OccupancyStudy
BIN_EXEC0    = $(addprefix $(BIN)/,$(EXEC0) )

default: all

# compile sources
$(LIB)/%.o: %.C
	@echo .
	@echo ... compiling source: $< to $@
	$(CC) $(CCFLAGS) $< $(ROOTINCS) $(INC) -o $@

# 0
$(BIN_EXEC0):
	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC0): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS0) ) )


clean:
	@echo ... cleaning
	rm -f $(LIB)/*.o
	rm -f $(LIB)/*.so
	rm -f $(BIN_EXEC0)

$(LIB):
	mkdir -p $(LIB)

$(BIN):
	mkdir -p $(BIN)

installdirs: $(LIB) $(BIN)

selector: installdirs $(BIN_EXEC0)

all: installdirs $(BIN_EXEC0) 
