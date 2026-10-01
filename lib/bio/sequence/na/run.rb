# frozen_string_literal: true

#
# = bio/sequence/na/run.rb - tandem repeats (runs) of nucleic acid sequences
#
# Copyright::   Copyright (C) 2026 BioRuby Project
# License::     The Ruby License
#

require 'bio/sequence' unless defined?(Bio::Sequence)
require 'bio/sequence/na' unless defined?(Bio::Sequence::NA)

module Bio
  class Sequence
    class NA < String

      # A tandem repeat run.
      #
      # A run is a maximal interval of a nucleic acid sequence whose
      # substring has a smallest period +p+ with exponent at least two,
      # i.e. it contains at least two consecutive copies (possibly a
      # partial one) of its period.  Extending the interval by one base
      # to the left or to the right destroys periodicity.
      #
      # Coordinates are half-open: +start_position+ is included and
      # +end_position+ is excluded.  Runs may overlap one another.
      class Run

        # The smallest period (in bases) of the run.  The interval does
        # not admit a smaller period.
        attr_reader :period

        # Position of the first base of the run (zero-based, inclusive).
        attr_reader :start_position

        # Position immediately after the last base of the run
        # (zero-based, exclusive).
        attr_reader :end_position

        # Creates a new Run instance.  Usually not called directly;
        # use Bio::Sequence::NA#runs instead.
        # ---
        # *Arguments*:
        # * (required) _start_position_: Integer (zero-based, inclusive)
        # * (required) _end_position_: Integer (zero-based, exclusive)
        # * (required) _period_: Integer (smallest period, in bases)
        # *Returns*:: Run object
        def initialize(start_position, end_position, period)
          @start_position = start_position
          @end_position = end_position
          @period = period
        end

        # Length of the run in bases.
        # ---
        # *Returns*:: Integer
        def length
          @end_position - @start_position
        end
        alias size length

        # Exponent (number of period repetitions, possibly fractional)
        # of the run.  Always at least 2.0.
        # ---
        # *Returns*:: Float
        def exponent
          length.to_f / @period
        end

        # Converts the run to an array of
        # +[start_position, end_position, period]+.
        # ---
        # *Returns*:: Array
        def to_a
          [@start_position, @end_position, @period]
        end

        # Returns +true+ if +other+ represents the same interval with the
        # same period.
        def ==(other)
          other.is_a?(Run) && to_a == other.to_a
        end
        alias eql? ==

        # Hash value consistent with #==.
        def hash
          to_a.hash
        end

        # Human-readable representation of the run.
        # ---
        # *Returns*:: String
        def inspect
          "#<#{self.class}: #{to_a.inspect}>"
        end

        # IUPAC ambiguity codes are treated as repeat boundaries: a run
        # never crosses them.  'u' is accepted as an ordinary base.
        CONCRETE_BASE_BYTES = "atgcu".bytes.freeze

        # Sentinel byte separating the pattern from the text in the
        # Z-function arrays.  DNA and RNA only contain bytes above it.
        SEPARATOR = 0

        # Finds every maximal tandem repeat run of +sequence+ and returns
        # them sorted by half-open interval.  Runs may overlap.
        #
        # The sequence is split at IUPAC ambiguity codes (ambiguous
        # bases act as boundaries) and runs are searched independently
        # in each concrete segment with the Main-Lorentz divide and
        # conquer algorithm.  At every recursion level only a constant
        # number of Z-function arrays is computed; no pair of positions
        # or intervals is enumerated, so the running time is O(n log n).
        # ---
        # *Arguments*:
        # * (required) _sequence_: String (nucleic acid sequence)
        # *Returns*:: Array of Run objects
        def self.find_all(sequence)
          bytes = sequence.is_a?(String) ? sequence.bytes : sequence
          concrete = concrete_lookup
          results = []
          segment_start = nil
          index = 0
          while index < bytes.length
            if concrete[bytes[index]]
              segment_start = index if segment_start.nil?
            elsif segment_start
              append_segment(bytes, segment_start, index, results)
              segment_start = nil
            end
            index += 1
          end
          if segment_start
            append_segment(bytes, segment_start, bytes.length, results)
          end
          results.sort_by!(&:to_a)
          results
        end

        class << self
          private

          # Lookup table marking concrete (non-ambiguous) base bytes.
          def concrete_lookup
            @concrete_lookup ||= begin
              table = Array.new(256, false)
              CONCRETE_BASE_BYTES.each { |byte| table[byte] = true }
              table
            end
          end

          # Searches one boundary-free segment and keeps only maximal,
          # minimal-period intervals.
          def append_segment(bytes, segment_start, segment_end, results)
            segment = bytes[segment_start...segment_end]
            candidates = find_candidates(segment)
            length = segment.length
            candidates.each do |(start, finish), period|
              next unless maximal?(segment, start, finish, period, length)
              interval = segment[start...finish]
              next unless minimal_period(interval) == period
              results << new(segment_start + start,
                            segment_start + finish, period)
            end
          end

          # Periodicity must not survive a one-base extension on either
          # side of the segment.
          def maximal?(segment, start, finish, period, length)
            if start > 0 &&
               segment[start - 1] == segment[start - 1 + period]
              return false
            end
            if finish < length && segment[finish - period] == segment[finish]
              return false
            end
            true
          end

          # Main-Lorentz recursion.  Candidates is a map from a
          # half-open [start, finish] interval to the smallest period
          # reported for it.
          def find_candidates(data)
            candidates = {}
            scan_ranges(data, 0, data.length, candidates)
            candidates
          end

          # Divides the data at the cut and reports candidate chains
          # crossing that cut before recursing into both halves.
          def scan_ranges(data, range_start, range_end, candidates)
            length = range_end - range_start
            return if length < 2

            split = length / 2
            cut = range_start + split
            left = data[range_start...cut]
            right = data[cut...range_end]
            right_length = length - split
            context = chain_context(left, right, split, right_length)
            max_period = split < right_length ? split : right_length
            1.upto(max_period) do |period|
              add_chain_candidates(context, cut, period, candidates)
            end
            scan_ranges(data, range_start, cut, candidates)
            scan_ranges(data, cut, range_end, candidates)
          end

          # Precomputes the four Z-function arrays that measure longest
          # common prefixes of periodic chains crossing a cut.
          def chain_context(left, right, split, right_length)
            full = left + right
            reverse = full.reverse
            reverse_left = left.reverse
            z_b = if right_length >= 2
                    z_function(full[(split + 1)..-1] + [SEPARATOR] + full)
                  end
            {
              length: full.length,
              split: split,
              right_length: right_length,
              z_a: z_function(reverse[(right_length - 1)..-1] +
                               [SEPARATOR] + reverse),
              z_a2: z_function(reverse_left + [SEPARATOR] + reverse),
              z_b: z_b,
              z_b2: z_function(right + [SEPARATOR] + full)
            }
          end

          # Reports the (at most two) maximal matched-pair chains of the
          # given period crossing the current cut.
          def add_chain_candidates(context, cut, period, candidates)
            a = left_extension(context, period)
            a2 = context[:z_a2][context[:length] + 1 - period]
            b = right_extension(context, period)
            b2 = far_right_extension(context, period)

            # Chain passing through the pair at the cut: pairs
            # [1-a, b] in the j parameterisation.
            if a >= 1 && a + b >= period
              add_candidate(candidates, cut + 1 - a - period,
                            cut + b + 1, period)
            end
            # Chain crossing the cut but avoiding the pair at the cut:
            # pairs [period-a2, period-1+b2].
            return unless a2 >= 1 && b2 >= 1 && a2 <= period - 1 &&
                          a2 + b2 >= period
            add_candidate(candidates, cut - a2,
                          cut + period + b2, period)
          end

          # Number of matched pairs extending left from the cut pair.
          def left_extension(context, period)
            cap = context[:split] - period + 1
            value = context[:z_a][context[:length] + period + 1]
            [value, cap].min
          end

          # Number of matched pairs extending right from the cut pair.
          def right_extension(context, period)
            return 0 if context[:right_length] < 2
            cap = context[:right_length] - 1
            value = context[:z_b][context[:length] - period + 1]
            [value, cap].min
          end

          # Number of matched pairs extending right one period away.
          def far_right_extension(context, period)
            return 0 if period > context[:right_length] - 1
            cap = context[:right_length] - period
            value = context[:z_b2][context[:length] + period + 1]
            [value, cap].min
          end

          # Records an interval, keeping the smallest period seen.
          def add_candidate(candidates, start, finish, period)
            return if finish - start < 2 * period
            key = [start, finish]
            old = candidates[key]
            candidates[key] = period if old.nil? || period < old
          end

          # Smallest period of an array using one Z-function pass.
          def minimal_period(data)
            z = z_function(data)
            length = data.length
            period = 1
            while period < length
              return period if z[period] >= length - period
              period += 1
            end
            length
          end

          # Z-function: z[i] is the length of the longest common prefix
          # of +data+ and its suffix beginning at position i.
          def z_function(data)
            length = data.length
            z = Array.new(length, 0)
            window_left = 0
            window_right = 0
            index = 1
            while index < length
              value = z_initial_value(z, index, window_left, window_right)
              value += 1 while index + value < length &&
                              data[value] == data[index + value]
              z[index] = value
              if index + value > window_right
                window_left = index
                window_right = index + value
              end
              index += 1
            end
            z
          end

          # Bootstrapping value of the Z-function from its rightmost
          # match window.
          def z_initial_value(z, index, window_left, window_right)
            return 0 if index >= window_right
            distance = window_right - index
            previous = z[index - window_left]
            [distance, previous].min
          end
        end
      end # Run

    end # NA
  end # Sequence
end # Bio
