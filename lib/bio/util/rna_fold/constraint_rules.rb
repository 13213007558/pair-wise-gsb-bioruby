# frozen_string_literal: true

# bio/util/rna_fold/constraint_rules.rb - Constraint installation helpers
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    class Constraints
      # Private helpers that install and cross-validate pair constraints.
      module Rules
        private

        def install(forced_pairs, forbidden_pairs, unpaired)
          install_forced_pairs(forced_pairs)
          install_forbidden_pairs(forbidden_pairs)
          install_unpaired(unpaired)
        end

        def install_forced_pairs(pairs)
          pairs.each { |pair| install_forced_pair(normalized_pair(pair)) }
        end

        def install_forced_pair(pair)
          left, right = pair
          if @unpaired.key?(left) || @unpaired.key?(right)
            raise ConflictError, 'forced pair contains a forced-unpaired base'
          end
          unless add_forced_partner(left, right) &&
                 add_forced_partner(right, left)
            raise ConflictError, 'base assigned to more than one forced pair'
          end
          return if compatible_with_forced_pairs?(left, right)
          raise ConflictError, 'forced pairs cross each other'
        end

        def add_forced_partner(base, partner)
          existing = @partner[base]
          return existing == partner if existing
          @partner[base] = partner
          true
        end

        def compatible_with_forced_pairs?(left, right)
          @partner.each_with_index do |partner, index|
            next unless partner && index < partner
            next if left == index && right == partner
            return false if crosses?(left, right, index, partner)
          end
          true
        end

        def crosses?(left_a, right_a, left_b, right_b)
          (left_a < left_b && left_b < right_a && right_a < right_b) ||
            (left_b < left_a && left_a < right_b && right_b < right_a)
        end

        def install_forbidden_pairs(pairs)
          pairs.each { |pair| install_forbidden_pair(normalized_pair(pair)) }
        end

        def install_forbidden_pair(pair)
          left, right = pair
          if @partner[left] == right && @partner[right] == left
            raise ConflictError, 'pair is both forced and forbidden'
          end
          @prohibited[pair_key(left, right)] = true
        end

        def install_unpaired(positions)
          positions.each do |base|
            index = base.to_int
            validate_index(index)
            install_unpaired_index(index)
          end
        end

        def install_unpaired_index(index)
          if @partner[index]
            raise ConflictError,
                  'base is both forced-paired and forced-unpaired'
          end
          @unpaired[index] = true
        end

        def normalized_pair(pair)
          unless pair.respond_to?(:to_a) && pair.to_a.length == 2
            raise ArgumentError, "pair expected, got #{pair.inspect}"
          end
          left, right = pair.to_a.map(&:to_int).sort
          validate_index(left)
          validate_index(right)
          if left == right
            raise ArgumentError, 'a base cannot pair with itself'
          end
          [left, right]
        end

        def validate_index(index)
          return if index >= 0 && index < @length
          raise IndexError,
                "base index #{index} outside sequence of length #{@length}"
        end

        def pair_key(left, right)
          [left, right].freeze
        end
      end
    end
  end
end
