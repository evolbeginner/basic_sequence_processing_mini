#! /bin/env ruby


require 'getoptlong'

require 'bio'


#########################################################
tree_file = nil
bootstrap_min = 0


#########################################################
opts = GetoptLong.new(
  ['-i', '-t', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-(i|t)$/
      tree_file = value
    when /^-b$/
      bootstrap_min = value.to_f
  end
end


#########################################################
treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
tree = treeio.next_entry.tree
tree.options[:bootstrap_style] = :traditional

tree.each_node do |node|
  if node.respond_to?(:bootstrap) and not node.bootstrap.nil?
    if node.bootstrap < bootstrap_min
      node.bootstrap = nil
    end
  end
end

output = tree.output_newick
output.gsub!(/[\n\s]/, '')

puts output


