CXX=g++
CXXFLAGS=-Wall -Werror -Wextra -I../overlap

multiple-alignment-test: multiple_alignment.o multiple_alignment_test.o multiple_alignment_builder.o ../overlap/overlapper.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -g -o $@ $^

.o: %.cpp
	$(CXX) $(CXXFLAGS) -g -c $< -o $@

clean:
	rm *.o 
