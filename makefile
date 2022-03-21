CC = g++
CC_INCLUDE = -I/mnt/c/Users/kyoko/Desktop/LSDGM_2d



TARGETS = Main Main2 Main3 Main4 Main5 Main6 Main7 MainTest MainIso

all: $(TARGETS)

Main: main.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main2: main2.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main3: main3.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main4: main4.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main5: main5.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main6: main6.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main7: main7.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

Main8: Main8.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic

MainTest: mainTest.cpp $(wildcard *.h)
	$(CC) $< -Wl,-rpath=/opt/gcc/4.7.0/lib64 -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-variadic-macros

MainIso: MainIso.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic

MainIsoBox: MainIsoBox.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic

MainBiaxBox: MainBiaxBox.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic

MainShear: MainShear.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic

clean:
	rm -f $(TARGETS)
again: clean $(TARGETS)
