require 'open-uri'
local_fname = 'frankengenome.fa'
if File.exist?(local_fname) == false
	url = 'http://www.cs.jhu.edu/~langmea/resources/frankengene1.fasta'
	data = open(url).read
	File.open(local_fname, "w"){|file| file.write(open(url).read)}
end

sequence = String.new

#enter reads into array
File.open(local_fname, "r").each do |line|
	if line.include?('>') == false 
  		sequence = sequence + line.gsub(/\s+/, "")
  	end
end

aProbability = Hash.new(0)
aProbability['AT'] = sequence.scan(/AT/).count*1.0
aProbability['AA'] = sequence.scan(/AA/).count*1.0
aProbability['AC'] = sequence.scan(/AC/).count*1.0
aProbability['AG'] = sequence.scan(/AG/).count*1.0

number = aProbability['AT'] / (aProbability['AG'] + aProbability['AA'] + aProbability['AC'] + aProbability['AT'])
p number