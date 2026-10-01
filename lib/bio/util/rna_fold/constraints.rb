# frozen_string_literal: true

# bio/util/rna_fold/constraints.rb - Constraints for RNA folding
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require 'bio/util/rna_fold/constraint_rules'

module Bio
  class RNAFold
    # Stores and validates forced and forbidden base-pair constraints.
    class Constraints
      include Rules

      # Base-pair constraints conflict with each other or the sequence.
      class ConflictError < ArgumentError; end

      # @param length [Integer] sequence length
      # @param forced_pairs [Array<Array<Integer>>] pairs that must be present
      # @param forbidden_pairs [Array<Array<Integer>>] pairs that must be absent
      # @param unpaired [Array<Integer>] bases that must remain unpaired
      def initialize(length, forced_pairs: nil, forbidden_pairs: nil,
                     unpaired: nil)
        @length = length
        @partner = Array.new(length)
        @prohibited = {}
        @unpaired = {}
        install(Array(forced_pairs), Array(forbidden_pairs), Array(unpaired))
        freeze
      end

      # @return [Array<Integer>] frozen copy of the forced pairing partners
      def forced_partners
        @partner.dup.freeze
      end

      # @return [Array<Array<Integer>>] frozen forced pairs
      def forced_pairs
        pairs = []
        @partner.each_with_index do |partner, index|
          pairs << [index, partner] if partner && index < partner
        end
        pairs.freeze
      end

      # @return [Array<Array<Integer>>] frozen forbidden pairs
      def forbidden_pairs
        @prohibited.keys.sort.freeze
      end

      # @return [Array<Integer>] frozen bases forced to be unpaired
      def unpaired_positions
        @unpaired.keys.sort.freeze
      end

      # @return [Integer, nil] forced partner of +base+
      def forced_partner(base)
        @partner[base]
      end

      # @return [Boolean] whether +base+ participates in a forced pair
      def forced_paired?(base)
        !@partner[base].nil?
      end

      # @return [Boolean] whether +base+ is forced to stay unpaired
      def forced_unpaired?(base)
        @unpaired.key?(base)
      end

      # @return [Boolean] whether pair (i, j) is compatible with constraints
      def pair_allowed?(left, right)
        partner = @partner[left]
        return false if partner && partner != right
        partner = @partner[right]
        return false if partner && partner != left
        return false if @unpaired.key?(left) || @unpaired.key?(right)
        !@prohibited.key?(pair_key(left, right))
      end
    end
  end
end
