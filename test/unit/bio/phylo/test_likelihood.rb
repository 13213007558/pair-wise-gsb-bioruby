#
# = test/unit/bio/phylo/test_likelihood.rb - unit test for Bio::Phylo::Likelihood
#
# License::     The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/phylo/likelihood'
require 'bio/alignment'

module Bio
  class TestPhyloLikelihood < Test::Unit::TestCase
    # Builds a Bio::Tree from a list of
    # [parent_name, child_name, length] edges.  Internal node names
    # may be nil.  Returns [tree, {name => node}].
    def build_tree(edges)
      tree = Bio::Tree.new
      nodes = {}
      edges.each do |parent_name, child_name, length|
        parent = nodes[parent_name] ||= Bio::Tree::Node.new(parent_name)
        child = nodes[child_name] ||= Bio::Tree::Node.new(child_name)
        tree.add_edge(parent, child, Bio::Tree::Edge.new(length))
      end
      [tree, nodes]
    end

    # Builds a caterpillar (pectinate) tree with n tips.
    def build_caterpillar(n, length)
      edges = []
      (n - 1).times do |i|
        internal = "i#{i}"
        edges << [internal, "t#{i}", length]
        edges << [internal, i + 1 < n - 1 ? "i#{i + 1}" : "t#{n - 1}", length]
      end
      build_tree(edges)
    end

    def assert_in_delta_array(expected, actual, delta = 1.0e-9)
      assert_equal(expected.size, actual.size)
      expected.zip(actual) do |e, a|
        assert_in_delta(e, a, delta)
      end
    end

    # --- the substitution model ------------------------------------

    def test_jc7
    end

    def test_jc69_transition_probabilities
      model = Bio::Phylo::JC69.new
      assert_equal([0.25, 0.5, 0.75, 1.0].map { |_| 0.25 },
                   model.state_frequencies)
      # P(0) is exactly the identity matrix (deterministic
      # transmission on zero-length branches).
      p0 = model.transition_ß
    end
  end
end
