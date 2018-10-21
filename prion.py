#!/usr/local/bin/python

### prion.py - An extensible framework for LGP
### Copyright (C) 2018 Christopher J. Crowe
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program. If not, see <https://www.gnu.org/licenses/>

import random

# This is usually a minimal set of functions for the generated program to use.
# Can't imagine many reasons to not assume they will be included
# unless you imagine a good reason to do so, you need these for an actual system
# To really function at all.

# The "machine" (or processor), at its core is a simulated stack processor. In certain cases,
# stack-based machines demonstrate superior performance. In most of the remaineder, they
# seem to function no better or worse. Here is the naming scheme that will be utilized
# by the registers of this stack-based VM. The final implementation of the processor is uknown,
# so at this point in the design, we need to be careful about implementation decisions/constraints. 
# I want prion to almost boot-strap itself up so the priomordial functions are going to be
# instantiated as fully-realized and executable functions for the interpreter. I do want some
# flexibility at this point. But this is a very critical point in the base design so some
# aspects bear dwelling upon:
# -- At this point we can assume that we will have an array of stacks which will be utilized
#    as the registers of the processor
# -- These registers are capable of interpreting some intrinsic things or concepts. Such as the 
#    "push" and "pop" of the stack operations. 
# -- We also assume the inputs are a RO bus from which a single value is transmitted when
#    addressed by index

# Should you find any collisions in your local namespace for some odd reason, you can globally alter
# the naming conventions from here and the system can bootstrap itself out of you namespace
PRION_OUTPUT_NAME = "PRION_STACK"
PRION_INPUT_NAME  = "PRION_WIRE"
PRION_ERC_NAME    = "PRION_ERC"

# These next two will "feel" a tad strange, but the are part of the bootstapping process for this system
# The code will bootstrap itself everytime it comes up so please pay attention if you're changing anything
PRION_TAG_DELIMITER = "&"
PRION_OUTPUT_LBL    = PRION_TAG_DELIMITER + "OUT"
PRION_INPUT_LBL     = PRION_TAG_DELIMITER + "WIRE"
PRION_ERC_LBL       = PRION_TAG_DELIMITER + "ERC"
PRION_FUNCTION_LBL  = PRION_TAG_DELIMITER + "FUNC"
PRION_LABEL_ARRAY   = [PRION_OUTPUT_LBL, PRION_INPUT_LBL, PRION_FUNCTION_LBL]

# Now we actually start writing the code for our primordial functions

PRION_OPCODE_INT = -1
PRION_OPCODE_NOP = 0

# prion_INP creation here
PRION_OPCODE_INP   = PRION_OPCODE_NOP + 1
PRION_FUNCNAME_INP = "PRION_INP"
PRION_INP_FUNCTION_TEMP = """
def &FUNC(&OUT , &WIRE):
	&OUT.append(&WIRE)
	return
"""
PRION_INP_FUNCTION_DEF = PRION_INP_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_INP)
PRION_INP_FUNCTION_DEF = PRION_INP_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)
PRION_INP_FUNCTION_DEF = PRION_INP_FUNCTION_DEF.replace(PRION_INPUT_LBL, PRION_INPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_INP_FUNCTION_DEF)  #
###############################

# prion_PUSH creation here
PRION_OPCODE_PUSH        = PRION_OPCODE_INP + 1
PRION_FUNCNAME_PUSH      = "PRION_PUSH"
PRION_PUSH_FUNCTION_TEMP = """
def &FUNC(&OUT , &ERC):
	&OUT.append(&ERC)
	return
"""
PRION_PUSH_FUNCTION_DEF = PRION_PUSH_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_PUSH)
PRION_PUSH_FUNCTION_DEF = PRION_PUSH_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)
PRION_PUSH_FUNCTION_DEF = PRION_PUSH_FUNCTION_DEF.replace(PRION_ERC_LBL, PRION_ERC_NAME)


###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_PUSH_FUNCTION_DEF) #
###############################

# prion_POP creation here
PRION_OPCODE_POP        = PRION_OPCODE_PUSH + 1
PRION_FUNCNAME_POP      = "PRION_POP"
PRION_POP_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 0:
		&OUT.pop()
	if len(&OUT) == 0:
		&OUT.append(0)
	return
"""
PRION_POP_FUNCTION_DEF = PRION_POP_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_POP)
PRION_POP_FUNCTION_DEF = PRION_POP_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_POP_FUNCTION_DEF)  #
###############################


# Arithmetic Functions
# Now we build our more complex arithmetic &nbsp;prion functions<br />
# -- Addition
# -- Subtraction
# -- Multiplication
# -- Division
# These functions need to validate the stack and store and in some cases
# temporary data. Some, like division, will require additional care
# it must be protected from div/0 errors. Also we must ensure that
# all functions leave the stack they operated on in a consistent
# and meaningful state for the next operation (and they will not
# know which operation is theoretically coming next in their raw form)

# prion_ADD creation here
PRION_OPCODE_ADD        = PRION_OPCODE_POP + 1
PRION_FUNCNAME_ADD      = "PRION_ADD"
PRION_ADD_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 1:
		&OUT.append(&OUT.pop() + &OUT.pop())
		return
	if len(&OUT) == 0:
		&OUT.append(0)
		return
	return
"""
PRION_ADD_FUNCTION_DEF = PRION_ADD_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_ADD)
PRION_ADD_FUNCTION_DEF = PRION_ADD_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_ADD_FUNCTION_DEF)  #
###############################

# prion_SUB creation here
PRION_OPCODE_SUB        = PRION_OPCODE_ADD + 1
PRION_FUNCNAME_SUB      = "PRION_SUB"
PRION_SUB_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 1:
		&OUT.append(&OUT.pop() - &OUT.pop())
		return
	if len(&OUT) == 0:
		&OUT.append(0)
		return
	return
"""

PRION_SUB_FUNCTION_DEF = PRION_SUB_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_SUB)
PRION_SUB_FUNCTION_DEF = PRION_SUB_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_SUB_FUNCTION_DEF)  #
###############################

# prion_MULT creation here
PRION_OPCODE_MULT        = PRION_OPCODE_SUB + 1
PRION_FUNCNAME_MULT      = "PRION_MULT"
PRION_MULT_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 1:
		&OUT.append(&OUT.pop() * &OUT.pop())
		return
	if len(&OUT) == 0:
		&OUT.append(0)
		return
	return
"""

PRION_MULT_FUNCTION_DEF = PRION_MULT_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_MULT)
PRION_MULT_FUNCTION_DEF = PRION_MULT_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_MULT_FUNCTION_DEF) #
###############################

# prion_DIV creation here
PRION_OPCODE_DIV        = PRION_OPCODE_MULT + 1
PRION_FUNCNAME_DIV      = "PRION_DIV"
PRION_DIV_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 1:
		x1 = &OUT.pop() * 1.0
		x2 = &OUT.pop() * 1.0
		if x2 == 0.0:
			&OUT.append(x1)
			&OUT.append(x2)
			return
		&OUT.append(x1 / x2)
		return
	if len(&OUT) == 0:
		&OUT.append(0)
		return
	return
"""

PRION_DIV_FUNCTION_DEF = PRION_DIV_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_DIV)
PRION_DIV_FUNCTION_DEF = PRION_DIV_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_DIV_FUNCTION_DEF)  #
###############################

# prion_IF creation here
PRION_OPCODE_IF        = PRION_OPCODE_DIV + 1
PRION_FUNCNAME_IF      = "PRION_IF"
PRION_IF_FUNCTION_TEMP = """
def &FUNC(&OUT):
	if len(&OUT) > 2:
		x1 = &OUT.pop()
		x2 = &OUT.pop()
		x3 = &OUT.pop()
		if x1 > x2:
			&OUT.append(x1)
			return
		&OUT.append(x3)
		return
	if len(&OUT) == 0:
		&OUT.append(0)
		return
	return
"""

PRION_IF_FUNCTION_DEF = PRION_IF_FUNCTION_TEMP.replace(PRION_FUNCTION_LBL, PRION_FUNCNAME_IF)
PRION_IF_FUNCTION_DEF = PRION_IF_FUNCTION_DEF.replace(PRION_OUTPUT_LBL, PRION_OUTPUT_NAME)

###############################
# Indentaton matters here     #
# please be careful if        #
# this is moved around        #
###############################
exec(PRION_IF_FUNCTION_DEF)   #
###############################

# To do: Add boolean
# To do: Add trig
# To do: Add hyperbolics

########################################################################
# You can define your own base functions in the area below
########################################################################

########################################################################
# Here is where we package up our functional primitives nicely for
# the consumption of classes that would have no interest in the above
########################################################################
PRION_OXF = { PRION_OPCODE_INP:PRION_INP, PRION_OPCODE_PUSH:PRION_PUSH,
              PRION_OPCODE_POP:PRION_POP, PRION_OPCODE_ADD:PRION_ADD,
              PRION_OPCODE_SUB:PRION_SUB, PRION_OPCODE_MULT:PRION_MULT,
              PRION_OPCODE_DIV:PRION_DIV, PRION_OPCODE_IF:PRION_IF }

PRION_OXK = { PRION_OPCODE_INP:"INP", PRION_OPCODE_PUSH:"PUSH",
              PRION_OPCODE_POP:"POP", PRION_OPCODE_ADD:"ADD",
              PRION_OPCODE_SUB:"SUB", PRION_OPCODE_MULT:"MULT",
              PRION_OPCODE_DIV:"DIV", PRION_OPCODE_IF:"IF" }

PRION_SXO = { "pop":PRION_OPCODE_POP, "+":PRION_OPCODE_ADD,
              "-":PRION_OPCODE_SUB, "*":PRION_OPCODE_MULT,
              "/":PRION_OPCODE_DIV, "?":PRION_OPCODE_IF }

PRION_FSET = [ PRION_SXO[key] for key in PRION_SXO ]

PRION_ERC_INT = 0
PRION_ERC_DBL = 1
PRION_MUTDEF  = 5
PRION_O_IDX   = 0
PRION_S_IDX   = 1
PRION_I_IDX   = 2
PRION_E_IDX   = 2

PRION_F_SLCT  = 1
PRION_I_SLCT  = 2
PRION_E_SLCT  = 3

########################################################################
# The genome class is going to be what you will ultimately be using
# pretty much constanly by every other object in the system. I am keeping
# it in herently simple:
# -- It has a type of ERC (Ephemoral Random Constant)
# -- It has an intrisic idea of bouns for whatever the ERC requested is
# -- It is capable of generating a single line of LGP code for our stack
#    processor (a codon) and you have at least three available types:
#    [INPUT, ERC, FUNCTION]
# You can ask for a specific class of codon, or simply call the basic
# codon() to get one of them randomly (with equal weight)
#
# A valid codon follows the three template types. each of the *codon()
# functions will give the folowing types of returns
# -- i_codon() = input_codon = [OPCODE, STACK_IDX, INPUT_IDX]
# -- e_codon() = erc_codon   = [OPCODE, STACK_IDX, ERC_VAL]
# -- f_codon() = func_codon  = [OPCODE, STACK_IDX]
#
# I append IDX to the stack and input identifiers since what we could be
# looking at an array of inputs but also, if solving for multiple outputs,
# we can have an arrsy of stack registers (the value left of the top of
# the stack register once the sequence terminates will denote the "answer").
########################################################################

class Genome:
	def __init__(self, wcount, scount, erc_low=0, erc_high=1, fset=PRION_FSET, erc_type=PRION_ERC_INT):
		self.wcount   = wcount   # Number of wires (Inputs)
		self.scount   = scount   # Number of registes (Outputs)
		self.fset     = fset     # Functions that this genome can offer (Defaults to PRION_FSET)
		self.erc_type = erc_type # Data type of ERC (Defaults to Integer)
		self.erc_low  = erc_low  # Lower bound of ERC (Defaults to 0)
		self.erc_high = erc_high # Higher bound of ERC (Defaults to 1)

	# -- e_codon() = erc_codon   = [OPCODE, STACK_IDX, ERC_VAL]
	def e_codon(self):
		if self.erc_type == PRION_ERC_INT:
			return [PRION_OPCODE_PUSH, random.randint(0, self.scount - 1), random.randint(self.erc_low, self.erc_high)]
		if self.erc_type == PRION_ERC_FLOAT:
			return [PRION_OPCODE_PUSH, random.randint(0, self.scount - 1), random.uniform(self.erc_low, self.erc_high)]
		return [PRION_OPCODE_NOP] 

	# -- i_codon() = input_codon = [OPCODE, STACK_IDX, INPUT_IDX]
	def i_codon(self):
		return [PRION_OPCODE_INP, random.randint(0, self.scount - 1), random.randint(0, self.wcount - 1)]

	# -- f_codon() = func_codon  = [OPCODE, STACK_IDX]
	def f_codon(self):
		if len(self.fset) > 0:
			return [random.choice(self.fset), random.randint(0, self.scount - 1)]
		return [PRION_OPCODE_NOP]
	
	# return one of the three base codons with equal probability
	def codon(self):
		x = random.choice([PRION_F_SLCT, PRION_I_SLCT, PRION_E_SLCT])
		if x == PRION_E_SLCT:
			return self.e_codon()
		if x == PRION_I_SLCT:
			return self.i_codon()
		if x == PRION_F_SLCT:
			return self.f_codon()
		return [PRION_OPCODE_NOP]

########################################################################
# The sequencer is where you would put all of your enhanced filters for
# the genrated programs. This one is faily strightforward. Keep
# requesting codons at randome until the disred length is reached
# You can inherit and overload the sequence() function in derived classes
# to suit your problem domain as necessary
########################################################################
PRION_WEIGHT_DEF      = [1, 1, 1]
PRION_WEIGHT_ERC_IDX  = 0
PRION_WEIGHT_INP_IDX  = 1
PRION_WEIGHT_FUNC_IDX = 2
PRION_ERC_REQ         = 0
PRION_INP_REQ         = 1
PRION_FUNC_REQ        = 2

class Sequencer:
        def __init__(self, genome, weights=PRION_WEIGHT_DEF):
		self.genome = genome
                if weights != PRION_WEIGHT_DEF:
                        self.wheel  =  [PRION_INP_REQ for x in range(weights[PRION_WEIGHT_INP_IDX])]
                        self.wheel  += [PRION_ERC_REQ for x in range(weights[PRION_WEIGHT_ERC_IDX])]
                        self.wheel  += [PRION_FUNC_REQ for x in range(weights[PRION_WEIGHT_FUNC_IDX])]
                else:
                        self.wheel = [PRION_INP_REQ, PRION_ERC_REQ, PRION_FUNC_REQ]

        def sequence(self, length):
                rseq = []
                while len(rseq) < length:
                        reqtype = random.choice(self.wheel)
                        if reqtype == PRION_INP_REQ:
                                rseq.append(self.genome.i_codon())
                        elif reqtype == PRION_ERC_REQ:
                                rseq.append(self.genome.e_codon())
                        elif reqtype == PRION_FUNC_REQ:
                                rseq.append(self.genome.f_codon())
                        else:
                                rseq.append([PRION_OPCODE_NOP])
                return rseq


########################################################################
# The decoder is another simple base class provided mainly as a
# convenience to the user / researcher. It enables a sequence of valid
# codons to be slightly more readable. This of it as a reverse assembler
# for the virtual stack processor.
########################################################################
class Decoder:
	def __init__(self, genome):
		self.genome = genome

	def decode(self, prog):
		readable = ""
		for x in prog:
			readable += str(PRION_OXK[x[PRION_O_IDX]]) + "\t" + str(x[1:]) + "\n"
		return readable

########################################################################
# The processor is the base implementation of the virtual stack
# processor on which the generated sequences will be run. It is assumed
# the processor can tak an array of inputs (the size of which is contained
# in the genome) and an array of stack registers (also specified in the
# genome. The base processor takes a prion program and runs the sequence
# against stacks which start with the initial state [0].
# This means that a stack register in this virtual processor will always
# have an initial value
########################################################################
class Processor:
	def __init__(self, genome):
		self.genome = genome

	def evaluate(self, prog, inputs):
		stacks = [[0] for x in range(self.genome.scount)]
		for line in prog:
			if line[PRION_O_IDX]   == PRION_OPCODE_INT:
				pass
			elif line[PRION_O_IDX] == PRION_OPCODE_NOP:
				pass
			elif line[PRION_O_IDX] == PRION_OPCODE_INP:
				stacks[line[PRION_S_IDX]].append(inputs[line[PRION_I_IDX]])
			elif line[PRION_O_IDX] == PRION_OPCODE_PUSH:
				stacks[line[PRION_S_IDX]].append(line[PRION_E_IDX])
			else:
				PRION_OXF[line[PRION_O_IDX]](stacks[line[PRION_S_IDX]])
		return [x[-1] for x in stacks]

########################################################################
# The SPLICER CLASS  is where we start seeing more of the biological
# aspect of the prion object model.
#
# The base splicer performs a simple crossover operation. It makes the
# basic assumption that two sequences can be of different lengths.
# To avoid sizing problems it will ensure that it picks a cutpoint that
# will not exceed the bounds of either sequence. This assumes that both
# sequences have been generated with compatable genomes.
#
# Desired bahaviour with sequences of varying lengths can be hard to
# imagine in all cases from the library perspective, so I am leaving it
# tunable by providing a few flags
# -- PRION_TRIM: Generate a sequence of the minimum of the 2 lengths
# -- PRION_BALANCE: Generate a sequence of the average of the 2 lengths
# -- PRION_GROW: Generate a sequence of the minimum of the 2 lengths
# By default, the base splicer will always trim. This can be modified
# by passing a different flag to the splice() function
########################################################################
PRION_TRIM    = 0
PRION_BALANCE = 1
PRION_GROW    = 2

class Splicer:
    def __init__(self, genome):
        self.genome = genome

    def splice(self, prog1, prog2, mode=PRION_TRIM):

        s1 = len(prog1)
        s2 = len(prog2)

        if s1 == s2:
            cutpoint = random.randint(0, s1 - 1)
            return prog1[0:cutpoint] + prog2[cutpoint:]

        # If we get here, we are dealing with sequences that are not
        # of the same length
        if mode == PRION_TRIM:
            if s1 < s2:
                cutpoint = random.randint(0, s1 - 1)
                return prog1[0:cutpoint] + prog2[cutpoint:s1 - 1]
            else:
                cutpoint = random.randint(0, s2 - 1)
                return prog2[0:cutpoint] + prog1[cutpoint:s2 - 1]
        if mode == PRION_GROW:
            if s1 < s2:
                cutpoint = random.randint(0, s1 - 1)
                return prog1[0:cutpoint] + prog2[cutpoint:]
            else:
                cutpoint = random.randint(0, s2 - 1)
                return prog2[0:cutpoint] + prog1[cutpoint:]
        if mode == PRION_BALANCE:
            if s1 < s2:
                cutpoint = random.randint(0, s1 - 1)
                return prog1[0:cutpoint] + prog2[cutpoint:s1 + random.randint(s1, s2) - 1]
            else:
                cutpoint = random.randint(0, s2 - 1)
                return prog2[0:cutpoint] + prog1[cutpoint:s2 + random.randint(s2, s1) - 1]


########################################################################
# The mutator also shows some of its roots in biological/evolutionary
# processes. It is given a genome and a percentage p. If the
# random.uniform(0.0, 100.0) returns less than or equal to p, the current
# codon is replaced by a random one from the available genome
########################################################################
class Mutator:
	def __init__(self, genome):
		self.genome = genome

	def mutate(self, seq, p=PRION_MUTDEF):
		idx = 0
		seqsz = len(seq)
                rseq = []
		for codon in seq:
			if random.randint(0, 100) <= p:
				rseq.append(self.genome.codon())
                        else:
                                rseq.append(codon)
		return rseq

########################################################################
# The fitness class dtermines the overall score of a generated sequence.
# By default, the class will simply use a euclidan distance formula
# (although can be extended to cover many more).
# We will make the assumption that the "fitter" a program is, the closer
# the value will approach 1.0 (although overall fitness will be a number
# between (0.0 - 1.0]
# Fitness will always be positive in this context, so if the Fitness class
# is unable to calculate it for whatever reason, it will return -1
########################################################################
PRION_FIT_EUCLID = 0

class Fitness:
	def fitness(self, param, ideal, ftype=PRION_FIT_EUCLID):
		if ftype == PRION_FIT_EUCLID:
			return self.f_euclidian(param, ideal)
		return -1

	def f_euclidian(self, param, ideal):
		if len(param) != len(ideal):
			return -1
		tmp = []
		idx = 0
		while idx < len(param):
			tmp.append((param[idx] - ideal[idx]) * (param[idx] - ideal[idx]) * 1.0)
			idx += 1
		return (1.0 / (1.0 + sum(tmp)))

########################################################################
# The population is a conatainer class for sequences. It is included as
# a helper (meaning you need not use it if you would rather manage
# populations of programs directly. It will utilize tournament selection
# for designating programs for "breeding" (crossover).
# the individual sequences are stored as a list of 4 elements with the
# fitness score stored at IDX0 and the source stored at IDX1.
# The geneology will be stored at IDX3 and the sequences unique ID for
# the current run will be stored at IDX4
########################################################################
PRION_POP_SCORE_IDX  = 0
PRION_POP_SOURCE_IDX = 1
PRION_POP_GENY_IDX   = 2
PRION_POP_ID_IDX     = 3
PRION_POP_PARAM_IDX  = 0
PRION_POP_IDEAL_IDX  = 1

PRION_POP_DEFSIZE    = 5000

PRION_CENSUS_HIGHLANDER = 0
PRION_CENSUS_HIDX       = 1
PRION_CENSUS_HSCORE     = 2
PRION_CENSUS_GEN        = 3
PRION_CENSUS_MAXPOP     = 4
PRION_CENSUS_CURPOP     = 5

class Population:
        def __init__(self, genome, bdata, maxpop=PRION_POP_DEFSIZE, fit=Fitness()):
                self.genome     = genome
                self.maxpop     = maxpop
                self.gencount   = 0
                self.fit        = fit
                self.ID         = 0
                self.bdata      = bdata
                self.proc       = Processor(genome)
                self.holding    = []
                self.highlander = -1
                self.hidx       = -1
                self.hscore     = -1
                self.splicer    = Splicer(self.genome)
                
        def add(self, Seq):
                if len(self.holding) >= self.maxpop:
                        return
                member = [[-1], [[]], [-1], [0]]
                # Add source
                member[PRION_POP_SOURCE_IDX] = Seq
                # Assign run ID
                member[PRION_POP_ID_IDX] = self.ID
                self.ID += 1
                # Score member
                member[PRION_POP_SCORE_IDX] = self.s_score(Seq)
                self.holding.append(member)
                if self.hscore < member[PRION_POP_SCORE_IDX]:
                        self.hscore     = member[PRION_POP_SCORE_IDX]
                        self.highlander = member[PRION_POP_ID_IDX]
                        self.hidx       = len(self.holding) - 1
                return
        
        def s_score(self, Seq):
                datac = len(self.bdata) 
                if datac == 0:
                        return -1

                rsum = 0.0
                for entry in self.bdata:
                        rsum += self.fit.fitness(self.proc.evaluate(Seq, entry[PRION_POP_PARAM_IDX]),
                                                 entry[PRION_POP_IDEAL_IDX])
                return (rsum / datac)

        def popdump(self):
                return self.holding

        def census(self):
                d = {}
                d[PRION_CENSUS_HIGHLANDER] = self.highlander
                d[PRION_CENSUS_HIDX]       = self.hidx
                d[PRION_CENSUS_HSCORE]     = self.hscore
                d[PRION_CENSUS_GEN]        = self.gencount
                d[PRION_CENSUS_MAXPOP]     = self.maxpop
                d[PRION_CENSUS_CURPOP]     = len(self.holding)
                return d

        def thunderdome(self, player1, player2):
                if len(player1) != 1:
                        return max(self.thunderdome(player1[:(len(player1)/2)], player1[len(player1)/2:]),
                                   self.thunderdome(player2[:(len(player2)/2)], player2[len(player2)/2:]))
                if self.holding[player1[0]][PRION_POP_SCORE_IDX] > self.holding[player2[0]][PRION_POP_SCORE_IDX]:
                        return player1[0]
                return player2[0]
                        
        def evolve(self, tsize=32):
                next_holding = []
                next_holding.append(self.holding[self.hidx])
                self.hidx = 0
                self.gencount += 1

                while len(next_holding) < self.maxpop:
                        poolcount = len(self.holding)
                        players = [random.randint(0, poolcount - 1) for x in range(tsize)]
                        winner1_idx = self.thunderdome(players[:(len(players)/2)], players[len(players)/2:])
                        players = [random.randint(0, poolcount - 1) for x in range(tsize)]
                        winner2_idx = self.thunderdome(players[:(tsize/2)], players[(tsize/2):])
                        offspring = self.splicer.splice(self.holding[winner1_idx][PRION_POP_SOURCE_IDX],
                                                        self.holding[winner2_idx][PRION_POP_SOURCE_IDX])
                
                        member = [[-1], [[]], [self.gencount], [winner1_idx, winner2_idx]]
                        # Add source
                        member[PRION_POP_SOURCE_IDX] = offspring
                        # Assign run ID
                        member[PRION_POP_ID_IDX] = self.ID
                        self.ID += 1
                        # Score member
                        member[PRION_POP_SCORE_IDX] = self.s_score(offspring)
                        next_holding.append(member)
                        if self.hscore < member[PRION_POP_SCORE_IDX]:
                                self.hscore     = member[PRION_POP_SCORE_IDX]
                                self.highlander = member[PRION_POP_ID_IDX]
                                self.hidx       = len(next_holding) - 1
                self.holding = next_holding
                return
        
def PRION_LIBTEST():
	print("Running test harness")
	print("Please set PRION_NOISY to True at the top of this file for full bootstrap details")
	print("Please set PRION_LIBTEST to True at the top of this file for full bootstrap tests")
	print("---------------------------------------------")
	print("Testing Genome Class")
	g = Genome(5, 2, erc_low=-10, erc_high=10)
	progs = [] 
	print("creating list of 20 random ops")
	progs.append([g.codon() for x in range(20)])
	print(str(progs[0]))
	print("creating list of 20 random ERC ops")
	progs.append([g.e_codon() for x in range(20)])
	print(str(progs[1]))
	print("creating list of 20 random input ops")
	progs.append([g.i_codon() for x in range(20)])
	print(str(progs[2]))
	print("creating list of 20 random function ops")
	progs.append([g.f_codon() for x in range(20)])
	print(str(progs[3]))
	print("---------------------------------------------")

	print("Testing Sequencer Class")
	s = Sequencer(g)
	print("creating list of 500 random programs of 256 codons in length")
	progs = [s.sequence(256) for x in range(500)]
	print("---------------------------------------------")

	print("Testing Decoder Class")
	d = Decoder(g)
	print("Decoding Random Program")
	print(d.decode(random.choice(progs)))
	print("---------------------------------------------")

	print("Testing Processor Class")
	p = Processor(g)
	print("Setting inputs to: [100, 20, 43, 91, 16]")
	i = [100, 20, 43, 91, 16]
	pcounter2 = 0
	for prog in progs:
		print("PROGRAM - " + str(pcounter2) + " (output)" + str(p.evaluate(prog, i)))
		pcounter2 += 1
	print("---------------------------------------------")

	print("Testing Splicer Class")
	s = Splicer(g)
	print("First program")
	rand_p1 = random.choice(progs)
	print("Return val: " + str(p.evaluate(rand_p1, i)))
	print("Second &nbsp;program")
	rand_p2 = random.choice(progs)
	print("Return val: " + str(p.evaluate(rand_p2, i)))
	spliced_p = s.splice(rand_p1, rand_p2)
	print("Spliced program")
	print("Return val: " + str(p.evaluate(spliced_p, i)))
	print("---------------------------------------------")

	print("Testing Mutator Class")
	m = Mutator(g)
	print("Mutating spliced program")
	print("Running 25 mutattion passes @ 30% (etremely high in most cases)")
	for x in range(25):
		print("Pass - " + str(x))
		spliced_p = m.mutate(spliced_p, p=30)
		print("Return val: " + str(p.evaluate(spliced_p, i)))
	print("---------------------------------------------")

	print("Testing Fitness Class")
	f = Fitness()
	print("Randomly selecting 10 programs - testing aginst ideal [1, 2]")
	ideal = [1, 2]
	for x in range(10):
		print("Pass - " + str(x))
		result = p.evaluate(random.choice(progs), i)
		print("Return val: " + str(result))
		print("Fitness: " + str(f.fitness(result, ideal)))
	print("Testing fitness with correct answer")
	print("Fitness: " + str(f.fitness(ideal, ideal)))
	print("---------------------------------------------")
	
	print("Tests complete")
	print("---------------------------------------------")


