# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 2,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/db/newick'

module Bio
  class TestTreeRobinsonFoulds < Test::Unit::TestCase
    def tree(newick)
      Bio::Newick.new(newick).tree
    end

    def labels(splits)
      splits.map(&:labels)
    end

    def test_unrooted_distance_and_one_sided_splits
      tree1 = tree('((A,B),(C,D));')
      tree2 = tree('((A,C),(B,D));')
      result = tree1.robinson_foulds(tree2)

      assert_equal(2, result.distance)
      assert_equal([%w[A B]], labels(result.first_tree_splits))
      assert_equal([%w[A C]], labels(result.second_tree_splits))
    end

    def test_restricts_trees_to_common_leaves_and_compresses_paths
      tree1 = tree('((A,X,B):0.1,Y,(C,D):0.3);')
      tree2 = tree('((A,C):1,(B,D):2,Z);')
      result = tree1.robinson_foulds(tree2)

      assert_equal(%w[A B C D], result.common_leaves)
      assert_equal(2, result.distance)
      assert_equal([%w[A B]], labels(result.first_tree_splits))
      assert_equal([%w[A C]], labels(result.second_tree_splits))
    end

    def test_comparison_does_not_modify_input_trees
      tree1 = tree('((A,X,B):0.1,Y,(C,D):0.3);')
      tree2 = tree('((A,C):1,(B,D):2,Z);')
      before = [tree1.number_of_nodes, tree1.number_of_edges,
                tree2.number_of_nodes, tree2.number_of_edges]

      tree1.robinson_foulds(tree2, branch_length: true)
      after = [tree1.number_of_nodes, tree1.number_of_edges,
               tree2.number_of_nodes, tree2.number_of_edges]

      assert_equal(before, after)
    end

    def test_unrooted_distance_is_invariant_to_display_root
      reference = tree('((A,B),(C,D));')

      assert_equal(0, reference.robinson_foulds(tree('(A,B,(C,D));')).distance)
      assert_equal(0,
                   reference.robinson_foulds(
                     tree('(D,(C,(B,A)));')
                   ).distance)
    end

    def test_rooted_distance_uses_clades
      tree1 = tree('(A,(B,(C,D)));')
      tree2 = tree('(B,(A,(C,D)));')
      unrooted = tree1.robinson_foulds(tree2)
      rooted = tree1.robinson_foulds(tree2, rooted: true)

      assert_equal(0, unrooted.distance)
      assert_equal(2, rooted.distance)
      assert_equal([%w[B C D]], labels(rooted.first_tree_splits))
      assert_equal([%w[A C D]], labels(rooted.second_tree_splits))
    end

    def test_same_tree_different_newick_representations_have_zero_distance
      reference = tree('((A,B),(C,D));')

      assert_equal(0,
                   reference.robinson_foulds(
                     tree('((D,C),(B,A));')
                   ).distance)
      assert_equal(0,
                   reference.robinson_foulds(
                     tree('(A,(((D,C),B)));')
                   ).distance)
    end

    def test_branch_length_differences_sum_compressed_path_lengths
      tree1 = tree('((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);')
      tree2 = tree('((A:0.2,X:9):0.7,' \
        '((B:1,Y:2):0.1,(C:0.4,D:0.5):0.2):0.5);')
      result = tree1.robinson_foulds(tree2, branch_length: true)

      assert_equal(0, result.distance)
      assert_in_delta(0.7, result.branch_length_distance, 1e-12)
      shared = result.branch_length_differences
      assert_equal(1, shared.size)
      split = shared.keys.find { |item| item.labels == %w[A B] }
      assert_not_nil(split)
      assert_in_delta(0.7, shared[split], 1e-12)
    end

    def test_duplicate_leaf_names_are_rejected
      tree1 = tree('((A,B),(C,D));')
      tree2 = tree('((A,A),(C,D));')

      assert_raise(ArgumentError) do
        tree1.robinson_foulds(tree2)
      end
    end

    def test_fewer_than_two_common_leaves_are_rejected
      tree1 = tree('((A,X),(Y,Z));')
      tree2 = tree('((A,B),(C,D));')

      assert_raise(ArgumentError) do
        tree1.robinson_foulds(tree2)
      end
    end

    def test_distance_convenience_methods
      tree1 = tree('((A,B):0.1,(C,D):0.3);')
      tree2 = tree('((A,C):1,(B,D):2);')

      assert_equal(2, tree1.robinson_foulds_distance(tree2))
      assert_in_delta(3.4,
                      tree1.robinson_foulds_distance(tree2, false, true),
                      1e-12)
    end
  end
end
