# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), '..', '..',
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/dcj'

module Bio
  class TestDCJ < Test::Unit::TestCase
    def test_inversion_distance_and_scenario
      source = [[1, 2, 3]]
      target = [[1, -2, 3]]

      assert_equal(1, DCJ.distance(source, target))
      assert_scenario(source, target)
    end

    def test_linear_fusion_and_fission
      assert_equal(1, DCJ.distance([[1], [2]], [[1, 2]]))
      assert_equal(1, DCJ.distance([[1, 2]], [[1], [2]]))
      assert_scenario([[1], [2]], [[1, 2]])
      assert_scenario([[1, 2]], [[1], [2]])
    end

    def test_mixed_linear_and_circular_genomes
      source = [[1, 2], { blocks: [3, 4], circular: true }]
      target = [[1, -2, 3, -4]]

      assert_equal(3, DCJ.distance(source, target))
      assert_scenario(source, target)
    end

    def test_circular_start_is_representation_invariant
      source = [{ blocks: [1, 2, 3], circular: true }]
      rotated = [{ blocks: [3, 1, 2], circular: true }]
      target = [[1, -2, 3]]

      assert_equal(DCJ.distance(source, target),
                   DCJ.distance(rotated, target))
      assert_equal(DCJ.scenario(source, target).length,
                   DCJ.scenario(rotated, target).length)
      assert_scenario(rotated, target)
    end

    def test_telomeres_are_not_ordinary_adjacencies
      genome = DCJ::Genome.new([[1, 2]])

      assert_equal(1, genome.adjacencies.length)
      assert_equal(2, genome.telomeres.length)
      assert_equal(1, DCJ.distance([[1, 2]], [[1], [2]]))
    end

    def test_rejects_missing_duplicate_and_invalid_blocks
      assert_raise(ArgumentError) { DCJ.distance([[1, 2]], [[1]]) }
      assert_raise(ArgumentError) { DCJ::Genome.new([[1, -1]]) }
      assert_raise(ArgumentError) { DCJ::Genome.new([[1, 0]]) }
      assert_raise(ArgumentError) { DCJ::Genome.new([[1, '2']]) }
      assert_raise(ArgumentError) { DCJ::Genome.new([[]]) }
      assert_raise(ArgumentError) { DCJ::Genome.new('invalid') }
      assert_raise(ArgumentError) do
        DCJ::Genome.new([{ blocks: [1], circular: 'yes' }])
      end
    end

    def test_exhaustive_state_graph_for_small_instances
      [2, 3].each do |count|
        states = all_matchings(endpoints(count))
        assert_equal(matching_count(2 * count), states.length)

        states.each do |source|
          states.each do |target|
            source_genome = genome_from_matching(source, count)
            target_genome = genome_from_matching(target, count)
            expected = bfs_distance(source, target, count)

            assert_equal(expected,
                         DCJ.distance(source_genome, target_genome))
            assert_equal(expected,
                         DCJ.scenario(source_genome, target_genome).length)
          end
        end
      end
    end

    private

    def assert_scenario(source, target)
      source_genome = DCJ::Genome.new(source)
      target_genome = DCJ::Genome.new(target)
      current = source_genome

      DCJ.scenario(source_genome, target_genome).each do |operation|
        current = current.apply(operation)
      end

      assert_equal(target_genome.adjacencies, current.adjacencies)
    end

    def endpoints(count)
      (1..count).flat_map do |id|
        [DCJ::Endpoint.new(id, :t), DCJ::Endpoint.new(id, :h)]
      end
    end

    def matching_count(vertex_count)
      (2..vertex_count).inject([1, 1]) do |counts, vertices|
        counts << counts[vertices - 1] + (vertices - 1) * counts[vertices - 2]
      end.last
    end

    def all_matchings(points)
      return [[]] if points.empty?

      first = points.first
      rest = points.drop(1)
      pairs = rest.each_with_index.flat_map do |point, index|
        remaining = rest.reject.with_index { |_item, i| i == index }
        all_matchings(remaining).map do |matching|
          [[first, point].sort] + matching
        end
      end
      (all_matchings(rest) + pairs).map(&:sort).uniq
    end

    def genome_from_matching(matching, count)
      DCJ::Genome.from_adjacencies(matching.to_set, (1..count).to_a)
    end

    def bfs_distance(source, target, count)
      queue = [source]
      distances = { source => 0 }

      until queue.empty?
        current = queue.shift
        return distances[current] if current == target

        legal_operations(current, count).each do |neighbor|
          next if distances.key?(neighbor)

          distances[neighbor] = distances[current] + 1
          queue << neighbor
        end
      end
    end

    def legal_operations(state, count)
      used = state.flatten
      free = endpoints(count) - used
      operations = state.map { |edge| state - [edge] }
      operations += free.combination(2).map do |pair|
        state + [pair.sort]
      end

      state.each do |edge|
        free.each do |point|
          edge.each do |endpoint|
            operations << (state - [edge]) + [[endpoint, point].sort]
          end
        end
      end

      state.combination(2) do |first, second|
        operations += reconnections(state, first, second)
      end

      operations.map(&:sort).uniq.reject { |operation| operation == state }
    end

    def reconnections(state, first, second)
      a, b = first
      c, d = second
      rest = state - [first, second]
      [
        rest + [[a, c].sort, [b, d].sort],
        rest + [[a, d].sort, [b, c].sort]
      ]
    end
  end
end
