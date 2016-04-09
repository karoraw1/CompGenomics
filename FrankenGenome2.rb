#!/usr/bin/env ruby

#these are the necessary files for training
require 'open-uri'
file1 = 'frankengenome.fa'
file2 = 'training.data'
file3 = 'test.data'

#go find them, and drop them locally for later use

if File.exist?(file1) == false
	url1 = 'http://www.cs.jhu.edu/~langmea/resources/frankengene1.fasta'
	File.open(file1, "w"){|file| file.write(open(url1).read)}
end

if File.exist?(file2) == false
	url2 = 'http://www.cs.jhu.edu/~langmea/resources/trainingData1.txt'
	File.open(file2, "w"){|file| file.write(open(url2).read)}
end

if File.exist?(file3) == false
	url3 = 'http://www.cs.jhu.edu/~langmea/resources/testData1.txt'
	File.open(file3, "w"){|file| file.write(open(url3).read)}
end

#import them into cache for right now

frankenSeq = String.new

File.open(file1, "r").each do |line|
	if line.include?('>') == false
  		frankenSeq << line.gsub(/\s+/, "")
  	end
end

trainingData = String.new
IO.foreach(file2) do |line|
	trainingData << line.gsub(/\s+/, "")
end

testData = String.new
IO.foreach(file3) do |line|
	testData << line.gsub(/\s+/, "")
end

#what would you have me do with them?

task = String.new
textToAnalyze = String.new
offset = 0

if ARGV[1].nil?
	textToAnalyze = frankenSeq
elsif ARGV[1] == "testData"
	offset = testData.length
	textToAnalyze = frankenSeq[offset..((offset+testData.length)-1)]
else
	textToAnalyze = frankenSeq
end

if ARGV[0].nil?
	task = "writeAnswer"
elsif ARGV[0] == "printResults" or ARGV[0] == "printScore" or ARGV[0] == "printSausage"
	task = ARGV[0]
else
	task = "writeAnswer"
end

#second order frequencies if necessary, however appropriate transition hash tables are not setup

#transitionDimers = Hash.new
#transitionDimers["AB"] = Hash.new(0)
#transitionDimers["BA"] = Hash.new(0)
#transitionDimers["BB"] = Hash.new(0)
#transitionDimers["AA"] = Hash.new(0)

#bin the training data by genome of origin into two strings
genomeA = String.new
genomeB = String.new

#keep track of transitions in order to split genomeX into array contig strings & to train transition matrix
checkpoints = [0]
transitionMatrix = Hash.new(0)

for i in 0..trainingData.length-1
	key = String.new

	#incorrect 2nd order transition matrix training

	#if i > 0
	#	if trainingData[i-1..i] == "00"
	#		key = "AA"
	#	elsif trainingData[i-1..i] == "11"
	#		key = "BB"
	#	elsif trainingData[i-1..i] == "01"
	#		key = "AB"
	#	elsif trainingData[i-1..i] == "10"
	#		key = "BA"
	#	end
	#	if transitionDimers[key].has_key?(frankenSeq[i])
	#		transitionDimers[key][frankenSeq[i]]+=1
	#	else
	#		transitionDimers[key][frankenSeq[i]]=1.0
	#	end
	#end

	#correct 1st order transition matrix training

	if i > 0
		thisState = trainingData[i]
		lastState = trainingData[i-1]
		if lastState != thisState
			checkpoints.push i
		end
		#counts transitions between states i.e. matrix A
		transitionMatrix[lastState+thisState] +=1
	end

	#assigns training data sequence to each genome
	if trainingData[i] == '0'
		genomeA << frankenSeq[i]
	end

	if trainingData[i] == '1'
		genomeB << frankenSeq[i]
	end
end

#translates transition matrix from counts into log probabilities in new hash
transitionMatrix2 = Hash.new

#denominators
zeros = genomeA.length*1.0
ones = genomeB.length*1.0

#also switch key nomenclature to letters
transitionMatrix.each do |key, value|
	if key == "00"
		transitionMatrix2["AA"] = Math.log2(value/zeros)
	end
	if key == "01"
		transitionMatrix2["AB"] = Math.log2(value/zeros)
	end
	if key == "10"
		transitionMatrix2["BA"] = Math.log2(value/ones)
	end
	if key == "11"
		transitionMatrix2["BB"] = Math.log2(value/ones)
	end

end

#unfinished second order transitions
#transitionDimers.each do |key, val|
#	denominator = val.values.inject(:+)
#	transitionDimers[key].each do |key2, val2|
#		transitionDimers[key][key2] = Math.log2(val2/denominator)
#	end
#end

#adds final checkpoint for calculating contig lengths
checkpoints.push trainingData.length

contigsA = []
contigsB = []

#turns each checkpoint index value into a contig length by subtracting subsequent indices and then assigning them to a genome
for i in 1..checkpoints.length-1
	contigLength = checkpoints[i]-checkpoints[i-1]
	if i.odd?
		contigsA << contigLength
	else
		contigsB << contigLength
	end
end

#cuts each genome from a continuous string, into the contigs observed in the frankengenome
contigAarray = Array.new
contigBarray = Array.new

contigsA.each do |contig|
	contigAarray.push genomeA.slice!(0..contig)
end

contigsB.each do |contig|
	contigBarray.push genomeB.slice!(0..contig)
end

#tallies mono, di, and trinucleotides observed in Genome A
#TODO: make this task a discrete function

onemerFreqA = Hash.new(0)
twomerFreqA = Hash.new(0)
threemerFreqA = Hash.new(0)

contigAarray.each do |contig|
	for i in 0..contig.length-1
		if onemerFreqA.has_key?(contig[i])
			onemerFreqA[contig[i]] += 1
		else
			onemerFreqA[contig[i]] = 1.0
		end
		if i < contig.length-2
			if twomerFreqA.has_key?(contig[i..i+1])
				twomerFreqA[contig[i..i+1]] += 1
			else
				twomerFreqA[contig[i..i+1]] = 1.0
			end
		end
		if i < contig.length-3
			if threemerFreqA.has_key?(contig[i..i+2])
				threemerFreqA[contig[i..i+2]] += 1
			else
				threemerFreqA[contig[i..i+2]] = 1.0
			end
		end
	end
end

#tallies mono, di, and trinucleotides observed in Genome B
#(only onemers used)

onemerFreqB = Hash.new(0)
twomerFreqB = Hash.new(0)
threemerFreqB = Hash.new(0)
contigBarray.each do |contig|
	for i in 0..contig.length-1
		if onemerFreqB.has_key?(contig[i])
			onemerFreqB[contig[i]] += 1
		else
			onemerFreqB[contig[i]] = 1.0
		end
		if i < contig.length-2
			if twomerFreqB.has_key?(contig[i..i+1])
				twomerFreqB[contig[i..i+1]] += 1
			else
				twomerFreqB[contig[i..i+1]] = 1.0
			end
		end
		if i < contig.length-3
			if threemerFreqB.has_key?(contig[i..i+2])
				threemerFreqB[contig[i..i+2]] += 1
			else
				threemerFreqB[contig[i..i+2]] = 1.0
			end
		end
	end
end

#turns trinucleotide counts into log probabilities for emission matrix, can also print them for debugging

=begin
threemersA = threemerFreqA.values.inject(:+)*1.0
threemersB = threemerFreqB.values.inject(:+)*1.0
threemerProbsA = Hash.new
threemerProbsB = Hash.new
threemerFreqA.each do |key, value|
	threemerProbsB[key] = Math.log2(threemerFreqB[key]/threemersB)
	threemerProbsA[key] = Math.log2(value/threemersA)
end

transitionDimers.each do |key, val|
	p key
	val.each do |key2, val2|
		print key2, " ", val2, " "
	end
end
threemerProbsB.each do |key, value|
	print key, ": "
	print value, " ", threemerProbsA[key], "\n"


#these are second order emission probabilities
#the "Loose" and "Strict" refer to whether the emissions are calculated using
#"Loose" = (all occurences of "CC" / all occurences of "CN")
#"Strict" = (all occurences of "CC" / all occurences of "C")
#they give very slighly different values, because Strict will always have a slightly larger denominator

genomeAemissions_loose = Hash.new
genomeAemissions_strict = Hash.new
genomeBemissions_loose = Hash.new
genomeBemissions_strict = Hash.new

twomerFreqA.each do |dimer, count|
	denominatorA = 0
	twomerFreqA.each do |key, count2|
		if key[0] == dimer[0]
			denominatorA = denominatorA + count2
		end
	end
	genomeAemissions_loose[dimer] = Math.log2((count*1.0)/denominatorA)
	genomeAemissions_strict[dimer] = Math.log2((count*1.0)/onemerFreqA[dimer[0]])
end

twomerFreqB.each do |dimer, count|
	denominatorB = 0
	twomerFreqB.each do |key, count2|
		if key[0] == dimer[0]
			denominatorB = denominatorB + count2
		end
	end
	genomeBemissions_loose[dimer] = Math.log2((count*1.0)/denominatorB)
	genomeBemissions_strict[dimer] = Math.log2((count*1.0)/onemerFreqB[dimer[0]])
end

#optional printing of 2nd order emissions
if task == "print2ndorder"
	print "           GenA_L\t             GenB_L\n"
	genomeAemissions_loose.each do |dimer, prob|
		print "  ", dimer[0], "|", dimer[1], ": ", prob, "\t",genomeBemissions_loose[dimer], "\n"
	end
	p transitionMatrix2
end
=end

#first order emissions for A
onemerProbsA = Hash.new
onemerFreqA.each do |base, count|
	onemerProbsA[base] = Math.log2((count*1.0)/(onemerFreqA.values.inject(:+)))
end
#first order emissions for B
onemerProbsB = Hash.new
onemerFreqB.each do |base, count|
	onemerProbsB[base] = Math.log2((count*1.0)/(onemerFreqB.values.inject(:+)))
end

#time to test the test data


#a path for each genome
viterbiPath = Hash.new
viterbiPath['A'] = [Math.log2(contigsA.inject(:+)/(trainingData.length-1.0))]
viterbiPath['B'] = [Math.log2(contigsB.inject(:+)/(trainingData.length-1.0))]

#keeps track of whether switching or staying was more likely for a given path
transitionCounter = Hash.new
transitionCounter["A"] = []
transitionCounter["B"] = []

for i in 1..textToAnalyze.length-1

	#emissionA = threemerProbsA[frankenSeq[offset+i-2..offset+i]]
	#emissionA = genomeAemissions_strict[frankenSeq[offset+i-1..offset+i]]
	#transitionAstay = transitionDimers['AA'][offset+i-1]+viterbiPath['A'][i-1]
	#transitionAswitch = transitionDimers['AB'][offset+i-1]+viterbiPath['B'][i-1]

	#add log probabilities instead of multiplying
	emissionA = onemerProbsA[frankenSeq[offset+i]]
	transitionAstay = transitionMatrix2['AA']+viterbiPath['A'][i-1]
	transitionAswitch = transitionMatrix2['AB']+viterbiPath['B'][i-1]

	if transitionAstay > transitionAswitch
		transitionCounter['A'][i-1] = 'stay'
	else
		transitionCounter['A'][i-1] = 'switch'
	end

	#emissionB = threemerProbsB[frankenSeq[offset+i-2..offset+i]]
	#emissionB = genomeBemissions_strict[frankenSeq[offset+i-1..offset+i]]
	#transitionBstay = transitionDimers['BB'][offset+i-1]+viterbiPath['B'][i-1]
	#transitionBswitch = transitionDimers['BA'][offset+i-1]+viterbiPath['A'][i-1]

	emissionB = onemerProbsB[frankenSeq[offset+i]]
	transitionBstay = transitionMatrix2['BB']+viterbiPath['B'][i-1]
	transitionBswitch = transitionMatrix2['BA']+viterbiPath['A'][i-1]

	if transitionBstay > transitionBswitch
		transitionCounter['B'][i-1] = 'stay'
	else
		transitionCounter['B'][i-1] = 'switch'
	end

	viterbiPath['A'][i] = emissionA+[transitionAstay, transitionAswitch].max
	viterbiPath['B'][i] = emissionB+[transitionBstay, transitionBswitch].max

	#debugging command
	if task == "printSausage"
		if i == 1
			puts ""
			print "char \t emission           \t stay           \t switch\n"
			print "         Paths  A: ", viterbiPath['A'][0], "\t", "B: ", viterbiPath['B'][0], "\n"
		end
		print frankenSeq[offset+i-1..offset+i], "\t", emissionA, "\t", transitionAstay, "\t", transitionAswitch, "\n"
		print " \t", emissionB, "\t", transitionBstay, "\t", transitionBswitch, "\n"
		print "         Paths  A: ", viterbiPath['A'][i], "\t", "B: ", viterbiPath['B'][i], "\n"
		puts ""
	end

end

#
backtrace = viterbiPath['A'].length-1
#initialize backtrace vector with first value
oneTruePath = String.new
if viterbiPath['A'][-1] > viterbiPath['B'][-1]
	oneTruePath << '0'
else
	oneTruePath << '1'
end

#backrace through path
for i in 1..backtrace
		#if we are in Genome A...
		if oneTruePath[i-1] == '0'
			#...and we stayed, predict genome A
			if transitionCounter["A"][-i] == 'stay'
				oneTruePath << '0'
			end
			#...if we switched, predict genome B
			if transitionCounter["A"][-i] == 'switch'
				oneTruePath << '1'
			end
		end
		#if we are instead in genome B...
		if oneTruePath[i-1] == '1'
			#...and we did not switch, predict more B
			if transitionCounter["B"][-i] == 'stay'
				oneTruePath << '1'
			end
			#...otherwise predict a switch to A
			if transitionCounter["B"][-i] == 'switch'
				oneTruePath << '0'
			end
		end
end

#since we backtraced, our parsimonious path is in reverse order
oneTruePath.reverse!

#what kind of answer shall we print?

if task == "printScore"
	truepositive = 0
	oneTruePath.each_char.each_with_index do |call, idx|
		if oneTruePath[idx] == testData[idx]
			truepositive += 1
		end
	end
	p (truepositive/50000.0)
end

if task == "printPath"
	print '   A                     B               path testData', "\n" 
	oneTruePath.each_char.each_with_index do |char, idx|
		print idx+1, ': ', viterbiPath['A'][idx], "  ", viterbiPath['B'][idx], "  ", char, '  ', testData[idx], "\n"
	end
end

if task == "printResults"
	for i in 0..(trainingData.length/100)-1
		print " ", oneTruePath[0+(i*100)..99+(i*100)],"\n"
		print " ", testData[0+(i*100)..99+(i*100)],"\n\n"
	end
end


if task == "writeAnswer"
	File.open("HMMResults.txt", "w") { |f|
		oneTruePath.each_char.each_with_index do |digit,index|
			f << digit
		end
	}
end


