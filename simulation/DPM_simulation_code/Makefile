CC = gcc
CFLAGS = -O2 -Wall -g

UTILDIR = utility

C_FILES = revisit_optimize_multidrug_responses_trimrates_strategies_index1.c mat_utility.c basic_utility.c drug_response_utility.c drug_response_utility2.c drug_response_utility3.c drug_response_utility4.c drug_response_utility5.c

INCLUDE_DIR = 
INCLUDE_FILES = db.h phylogeny.h phylogeny2.h
H_FILES = stdio.h stdlib.h string.h math.h
O_FILES = $(C_FILES:.c=.o)

EXECUTABLE = revisit_optimize_multidrug_responses_trimrates_strategies_index1

$(EXECUTABLE): $(O_FILES)
	$(CC) $(CFLAGS) -o $@ $(O_FILES) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm $(EXECUTABLE) $(O_FILES)
