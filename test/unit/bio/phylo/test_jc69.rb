
load Pathname.new(File.join(File.dirname(__FILE__), %w[.. .. .. bioruby_test_helper.rb])).cleanpath.to_s

require 'test/unit'
require 'bio'

module Bio
  class TestPhyloJC69 < Test::Unit::TestCase
    IUPAC_STATES = {
      'A' => [0], 'C' => [1], 'G' => [2], 'T' => [3], 'U' => [3],
      'M' => [0, 1], 'R' => [0, 2], 'W' => [0, 3], 'S' => [1, 2],
      'Y' => [1, 3], 'K' => [2, 3], 'V' => [0, 1, 2],
      'H' => [0, 1, 3], 'D' => [0, 2, 3], 'B' => [1, 2, 3],
      'N' => [0, 1, 2, 3]
    }.freeze

    def test_matches_enumeration_on_a_multifurcating_tree
      alignment = {
        'a' => 'ACGTNRMW',
        'b' => 'AGRYM-C?',
        'c' => 'ATWKSB-.',
        'd' => 'ACGVHDN-',
        'e' => 'AGTBMCWN'
      }
      tree = star_tree(alignment.keys, 0.37)
      categories = [[0.0, 0.2], [0.8, 0.3], [2.4, 0.5]]

      assert_enumerated_likelihood(tree, alignment, categories)
    end

    def test_matches_enumeration_on_a_binary_tree
      alignment = {
        'a' => 'ACGTN',
        'b' => 'AGRYM',
        'c' => 'ATWKS',
        'd' => 'ACGVH'
      }
      tree = binary_tree(0.13, 0.0, 0.27, 0.19, 0.33, 0.08)
      categories = [[0.25, 1.0], [1.5, 2.0]]

      assert_enumerated_likelihood(tree, alignment, categories)
    end

    def test_zero_length_branch_is_deterministic_for_equal_states
      tree = Bio::Tree.new
      a = make_node('a')
      b = make_node('b')
      tree.add_edge(a, b, Bio::Tree::Edge.new(0.0))

      assert_in_delta(-Math.log(4.0),
                      Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'b' => 'A'),
                      1e-12)
      assert_equal(-Float::INFINITY,
                   Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'b' => 'C'))
    end

    def test_single_tip_tree_uses_root_prior
      tree = Bio::Tree.new
      tree.add_node(make_node('a'))

      assert_in_delta(-Math.log(4.0),
                      Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A'),
                      1e-12)
      assert_equal(0.0,
                   Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'N'))
    end

    def test_all_missing_sites_have_unit_likelihood
      tree = star_tree(%w[a b c], 0.7)
      alignment = { 'a' => '?-', 'b' => '.N', 'c' => 'N?' }

      assert_equal([0.0, 0.0],
                   Bio::Phylo::JC69.site_log_likelihoods(tree, alignment))
      assert_equal(0.0, Bio::Phylo::JC69.log_likelihood(tree, alignment))
    end

    def test_weights_are_normalized
      tree = binary_tree(0.2, 0.3, 0.1, 0.5, 0.4, 0.0)
      alignment = { 'a' => 'ACG', 'b' => 'AGT', 'c' => 'ATC', 'd' => 'ACG' }
      unnormalized = [[1.0, 2.0], [3.0, 2.0]]
      normalized = [[1.0, 1.0], [3.0, 1.0]]

      assert_in_delta(Bio::Phylo::JC69.log_likelihood(tree, alignment, normalized),
                      Bio::Phylo::JC69.log_likelihood(tree, alignment, unnormalized),
                      1e-12)
    end

    def test_zero_weight_categories_are_ignored
      tree = star_tree(%w[a b c], 0.4)
      alignment = { 'a' => 'A', 'b' => 'C', 'c' => 'G' }
      with_zero_weight = [[1.0, 1.0], [5.0, 0.0]]

      assert_in_delta(Bio::Phylo::JC69.log_likelihood(tree, alignment),
                      Bio::Phylo::JC69.log_likelihood(tree, alignment,
                                                     with_zero_weight),
                      1e-12)
    end

    def test_long_tree_does_not_underflow
      labels = Array.new(1000) { |index| "leaf#{index}" }
      tree = star_tree(labels, 500.0)
      alignment = labels.each_with_object({}) { |label, hash| hash[label] = 'A' }
      result = Bio::Phylo::JC69.log_likelihood(tree, alignment)

      assert_predicate(result, :finite?)
      assert(result < Math.log(Float::MIN))
      assert_equal(0.0, 0.25 ** 1000)
      assert_in_delta(-1000.0 * Math.log(4.0), result, 1e-8)
    end

    def test_labels_must_match_alignment_rows
      tree = star_tree(%w[a b], 0.1)

      assert_raise(ArgumentError) do
        Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'c' => 'A')
      end
      assert_raise(ArgumentError) do
        Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'b' => 'A', 'c' => 'A')
      end
    end

    def test_accepts_bioruby_sequence_hash
      tree = star_tree(%w[a b c], 0.2)
      alignment = Bio::Alignment::SequenceHash[
        'a' => 'A', 'b' => 'C', 'c' => 'G'
      ]
      calculation = Bio::Phylo::JC69.calculate(tree, alignment)

      assert_kind_of(Bio::Phylo::JC69, calculation)
      assert_equal(1, calculation.site_log_likelihoods.length)
      assert_in_delta(calculation.site_log_likelihoods.first,
                      calculation.log_likelihood,
                      1e-15)
    end

    def test_rejects_invalid_branch_length
      tree = star_tree(%w[a b], nil)

      assert_raise(ArgumentError) do
        Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'b' => 'A')
      end
    end

    def test_rejects_unlabelled_alignment_and_invalid_nucleotide
      tree = star_tree(%w[a b], 0.1)

      assert_raise(ArgumentError) do
        Bio::Phylo::JC69.log_likelihood(tree,
                                        Bio::Alignment::SequenceArray['A', 'A'])
      end
      assert_raise(ArgumentError) do
        Bio::Phylo::JC69.log_likelihood(tree, 'a' => 'A', 'b' => 'X')
      end
    end

    private

    def make_node(name = nil)
      Bio::Tree::Node.new(name)
    end

    def add_edge(tree, parent, child, length)
      tree.add_edge(parent, child, Bio::Tree::Edge.new(length))
    end

    def star_tree(labels, length)
      tree = Bio::Tree.new
      root = make_node
      tree.root = root
      labels.each do |label|
        add_edge(tree, root, make_node(label), length)
      end
      tree
    end

    def binary_tree(root_a, root_b, a_length, b_length, c_length, d_length)
      tree = Bio::Tree.new
      root = make_node
      left = make_node
      right = make_node
      a = make_node('a')
      b = make_node('b')
      c = make_node('c')
      d = make_node('d')
      tree.root = root
      add_edge(tree, root, left, root_a)
      add_edge(tree, root, right, root_b)
      add_edge(tree, left, a, a_length)
      add_edge(tree, left, b, b_length)
      add_edge(tree, right, c, c_length)
      add_edge(tree, right, d, d_length)
      tree
    end

    def assert_enumerated_likelihood(tree, alignment, categories)
      expected = Array.new(alignment.values.first.length) do |site|
        weighted_likelihood = categories.sum do |rate, weight|
          weight * enumerated_site_likelihood(tree, alignment, site, rate)
        end
        weight_sum = categories.sum { |_rate, weight| weight }
        Math.log(weighted_likelihood / weight_sum)
      end
      actual = Bio::Phylo::JC69.site_log_likelihoods(tree, alignment, categories)

      expected.each_with_index do |value, site|
        assert_in_delta(value, actual[site], 1e-10, "site #{site}")
      end
      assert_in_delta(expected.sum, actual.sum, 1e-10)
    end

    def enumerated_site_likelihood(tree, alignment, site, rate)
      root = tree.nodes.find { |node| tree.out_degree(node) != 1 } || tree.nodes.first
      root_states = observed_states(tree, root, alignment, site)

      (0...4).sum do |state|
        next 0.0 unless root_states.nil? || root_states.include?(state)
        0.25 * enumerate_state(tree, alignment, site, rate, root, state, nil)
      end
    end

    def enumerate_state(tree, alignment, site, rate, node, state, parent)
      own_states = observed_states(tree, node, alignment, site)
      return 0.0 unless own_states.nil? || own_states.include?(state)

      tree.adjacent_nodes(node).reject { |neighbor| neighbor == parent }.reduce(1.0) do |result, child|
        edge = tree.get_edge(node, child)
        length = tree.get_edge_distance(edge)
        subtree = (0...4).sum do |child_state|
          child_states = observed_states(tree, child, alignment, site)
          next 0.0 unless child_states.nil? || child_states.include?(child_state)
          transition_probability(state, child_state, rate * length) *
            enumerate_state(tree, alignment, site, rate, child, child_state, node)
        end
        result * subtree
      end
    end

    def observed_states(tree, node, alignment, site)
      label = tree.get_node_name(node)
      return nil if label.nil?
      character = alignment.fetch(label, nil)&.slice(site)
      return nil if character.nil? || ['-', '.', '?'].include?(character.upcase)
      IUPAC_STATES.fetch(character.upcase)
    end

    def transition_probability(parent_state, child_state, scaled_length)
      survival = Math.exp(-4.0 * scaled_length / 3.0)
      if parent_state == child_state
        0.25 + 0.75 * survival
      else
        0.25 - 0.25 * survival
      end
    end
  end
end
