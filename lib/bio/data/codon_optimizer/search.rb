# frozen_string_literal: true

require 'bio/data/codon_optimizer/support'

module Bio
  class CodonOptimizer
    # Exact dynamic program for feasible codon assignments.
    # It examines every synonymous assignment through suffix/GC states, so
    # forbidden matches that cross codon boundaries are not repaired after a
    # greedy choice.
    class Search

      def initialize(choices, weights, forbidden_substrings, gc_range)
        @choices = choices
        @weights = weights
        @forbidden_substrings = forbidden_substrings
        @gc_range = gc_range
      end

      def run
        states = { ['', 0] => initial_state }
        @choices.each do |codons|
          states = StateTransition.new(
            states,
            codons,
            @weights,
            @forbidden_substrings,
            suffix_length
          ).next_states
          return [nil, forbidden_conflict] if states.empty?
        end
        StateSelector.new(states, @gc_range).best
      end

      private

      def initial_state
        { gc_count: 0, score: 0, codon: nil, parent: nil }
      end

      def suffix_length
        @suffix_length ||= @forbidden_substrings.empty? ? 0 :
          @forbidden_substrings.map(&:length).max - 1
      end

      def forbidden_conflict
        Support.conflict(
          'forbidden_substrings_impossible',
          'no synonymous coding avoids all forbidden substrings',
          forbidden_substrings: @forbidden_substrings
        )
      end

    end

    class StateTransition

      def initialize(states, codons, weights, forbidden_substrings,
                     suffix_length)
        @states = states
        @codons = codons
        @weights = weights
        @forbidden_substrings = forbidden_substrings
        @suffix_length = suffix_length
      end

      def next_states
        @codons.each_with_object({}) do |codon, next_states|
          @states.each do |(suffix, _gc), state|
            add_state(next_states, suffix, state, codon)
          end
        end
      end

      private

      def add_state(next_states, suffix, state, codon)
        combined = suffix + codon
        return if forbidden?(combined)

        key = [
          suffix_for(combined),
          state[:gc_count] + Support.gc_count(codon)
        ]
        score = state[:score] + @weights.fetch(codon, 0)
        replace(next_states, key, score, codon, state)
      end

      def replace(next_states, key, score, codon, state)
        existing = next_states[key]
        return if existing && existing[:score] >= score

        next_states[key] = {
          gc_count: key.last,
          score: score,
          codon: codon,
          parent: state
        }
      end

      def suffix_for(sequence)
        return '' if @suffix_length.zero?

        sequence[[sequence.length - @suffix_length, 0].max, @suffix_length]
      end

      def forbidden?(sequence)
        @forbidden_substrings.any? { |pattern| sequence.include?(pattern) }
      end

    end

    class StateSelector

      def initialize(states, gc_range)
        @states = states
        @gc_range = gc_range
      end

      def best
        candidates = feasible_states
        return [nil, gc_conflict] if candidates.empty?

        best = candidates.max_by { |state| state[:score] }
        [best.merge(codons: reconstruct(best)), nil]
      end

      private

      def feasible_states
        @states.values.select do |state|
          @gc_range.nil? || @gc_range.cover?(state[:gc_count])
        end
      end

      def reconstruct(state)
        codons = []
        while state && state[:codon]
          codons.unshift(state[:codon])
          state = state[:parent]
        end
        codons
      end

      def gc_conflict
        minimum, maximum = @states.keys.map(&:last).minmax
        Support.conflict(
          'gc_range_impossible',
          "achievable GC count is #{minimum}..#{maximum}, " \
            "but #{@gc_range.inspect} was requested",
          requested_min: @gc_range.begin,
          requested_max: @gc_range.end,
          achievable_min: minimum,
          achievable_max: maximum
        )
      end

    end

  end

end
