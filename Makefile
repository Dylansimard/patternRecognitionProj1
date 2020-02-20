TARGET = p1
LIBS = -lm #Math Library, just a placeholder
HEADERS = 
SRCS = as1.cpp boxmuller.cpp
OBJECTS := $(patsubst %.cpp,%.o,$(SRCS))
CXX = g++
CXX_FLAGS = -Wall -std=c++11

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f *.o 
	-rm -f $(TARGET)
	-rm -f *.txt
