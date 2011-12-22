CXX=g++
CXXFLAGS=-Wall -Werror -I../overlap

multiple-alignment-test: multiple_alignment.o multiple_alignment_test.o ../overlap/overlapper.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o 
