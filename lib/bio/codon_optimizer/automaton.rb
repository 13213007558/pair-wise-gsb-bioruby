# frozen_string_literal: true

# bio/codon_optimizer/automaton.rb - Aho-Corasick automaton for codon optimizer

module Bio
  class CodonOptimizer
    # Deterministic Aho--Corasick automaton over the DNA alphabet.
    class Automaton
      BASES = 'acgt'
      INDEX = BASES.each_char.with_index.to_h.freeze

      attr_reader :initial_state

      def initialize(patterns)
        @go = [[nil, nil, nil, nil]]
        @output = [[]]
        @failure = [0]
        patterns.each { |pattern| insert(pattern) }
        build
        @initial_state = 0
      end

      def transition(state, codon)
        codon.each_char do |base|
          state = @go[state][INDEX.fetch(base)]
          return nil unless @output[state].empty?
        end
        state
      end

      private

      def insert(pattern)
        state = 0
        pattern.each_char do |base|
          index = INDEX.fetch(base)
          unless @go[state][index]
            @go[state][index] = @go.length
            @go << [nil, nil, nil, nil]
            @output << []
          end
          state = @go[state][index]
        end
        @output[state] << pattern unless @output[state].include?(pattern)
      end

      def build
        queue = []
        BASES.each_char do |base|
          index = INDEX.fetch(base)
          child = @go[0][index]
          if child
            @failure[child] = 0
            queue << child
          else
            @go[0][index] = 0
          end
        end

        until queue.empty?
          state = queue.shift
          BASES.each_char do |base|
            index = INDEX.fetch(base)
            child = @go[state][index]
            if child
              queue << child
              fallback = @failure[state]
              until fallback.zero? || @go[fallback][index]
                fallback = @failure[fallback]
              end
              @failure[child] = @go[fallback][index] || 0
              @output[child] |= @output[@failure[child]]
            else
              @go[state][index] = @go[@failure[state]][index]
            end
          end
        end
      end
    end
  end
end
