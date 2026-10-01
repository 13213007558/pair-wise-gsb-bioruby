# frozen_string_literal: true

# bio/codon_optimizer/solver.rb - exact dynamic program

module Bio
  class CodonOptimizer
    # Exact dynamic program over automaton states and GC counts.
    class Solver
      Node = Struct.new(:score, :gc_count, :state, :codon, :predecessor,
                        keyword_init: true)

      def initialize(problem)
        @problem = problem
      end

      def run(codons, automaton, maximum_gc = problem.max_gc_count)
        layer = initial_layer(automaton)
        codons.each_with_index do |original, index|
          candidates = candidates(original, index, codons.length)
          layer = add_codon(layer, candidates, automaton, maximum_gc)
          break if layer.empty?
        end
        layer
      end

      def best(layer)
        layer.values.select { |node| node.gc_count >= problem.min_gc_count }
             .max_by(&:score)
      end

      def reconstruct(node)
        codons = []
        while node.predecessor
          codons.unshift(node.codon)
          node = node.predecessor
        end
        codons.join
      end

      private

      attr_reader :problem

      def initial_layer(automaton)
        root = Node.new(score: 0, gc_count: 0, state: automaton.initial_state,
                        codon: nil, predecessor: nil)
        { [root.state, root.gc_count] => root }
      end

      def candidates(original, index, codon_count)
        return [original] if problem.fixed?(original, index, codon_count)

        problem.table.revtrans(problem.table[original]).select do |codon|
          codon.match?(/\A[acgt]{3}\z/)
        end
      end

      def add_codon(layer, candidates, automaton, maximum_gc)
        candidates.each_with_object({}) do |codon, next_layer|
          weight = problem.weights.fetch(codon, 0)
          codon_gc = codon.count('gc')
          layer.each_value do |node|
            add_transition(node, codon, weight, codon_gc,
                           next_layer, automaton, maximum_gc)
          end
        end
      end

      def add_transition(node, codon, weight, codon_gc, next_layer,
                         automaton, maximum_gc)
        next_state = automaton.transition(node.state, codon)
        return unless next_state

        next_gc = node.gc_count + codon_gc
        return if maximum_gc && next_gc > maximum_gc

        next_score = node.score + weight
        key = [next_state, next_gc]
        return if next_layer[key] && next_layer[key].score >= next_score

        next_layer[key] = Node.new(score: next_score, gc_count: next_gc,
                                  state: next_state, codon: codon,
                                  predecessor: node)
      end
    end
  end
end
