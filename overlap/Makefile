CXX=g++
CXXFLAGS=-Wall -Werror -g

overlapper-test: overlapper.o overlapper-test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o 
