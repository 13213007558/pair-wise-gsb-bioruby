# frozen_string_literal: true

# = bio/util/rna_fold.rb - RNA partition function and pairing probabilities
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require 'bio/util/rna_fold/constants'
require 'bio/util/rna_fold/construction'
require 'bio/util/rna_fold/constraint_rules'
require 'bio/util/rna_fold/constraints'
require 'bio/util/rna_fold/energy_adapter'
require 'bio/util/rna_fold/partition_function'
require 'bio/util/rna_fold/validation'

module Bio
  # Pseudoknot-free RNA secondary structure partition function and
  # base-pair marginal (edge) probabilities.
  #
  # Unlike minimum-free-energy folding, the partition function weights every
  # non-crossing secondary structure with its Boltzmann factor and reports
  # how reliably each base pair or unpaired state is supported.
  #
  # Only additive pairing energies supplied by the caller are used:
  # +weight(pair) = exp(-beta * energy(i, j))+. No thermodynamic parameter
  # tables or external folding programs are required.
  #
  # == Example
  #
  #   energies = Hash.new { |hash, key| -1.0 }
  #   folder = Bio::RNAFold.new('GGGAAACCC', pairing_energy: energies)
  #   folder.log_partition_function
  #   folder.pair_probability(0, 8)
  #   folder.unpaired_probability(4)
  class RNAFold
    include Construction
    include PartitionFunction
    include Validation

    # Default minimum number of unpaired bases inside a hairpin loop.
    DEFAULT_MIN_HAIRPIN_LENGTH = 3

    # @return [String] sequence whose secondary structure is evaluated
    attr_reader :sequence

    # @return [Integer] minimum hairpin loop length
    attr_reader :min_hairpin_length

    # @return [Numeric] inverse Boltzmann factor applied to pairing energies
    attr_reader :beta

    # @return [Bio::RNAFold::Constraints] pairing constraints
    attr_reader :constraints

    # @param sequence [String, #to_s, #to_str] bases indexed from zero
    # @param pairing_energy [nil, #call, Hash, #[]] additive pair energies
    # @param min_hairpin_length [Integer] required unpaired bases per hairpin
    # @param beta [Numeric] inverse temperature, +1/(kT)+ in energy units
    # @param forced_pairs [Array<Array<Integer>>] pairs that must occur
    # @param forbidden_pairs [Array<Array<Integer>>] pairs that must not occur
    # @param unpaired [Array<Integer>] bases forced to remain unpaired
    def initialize(sequence, pairing_energy: nil,
                   min_hairpin_length: DEFAULT_MIN_HAIRPIN_LENGTH,
                   beta: 1.0, forced_pairs: nil, forbidden_pairs: nil,
                   unpaired: nil)
      @sequence = convert_sequence(sequence)
      @min_hairpin_length =
        validate_nonnegative_integer(min_hairpin_length, :min_hairpin_length)
      @beta = validate_positive_number(beta, :beta)
      @constraints = build_constraints(forced_pairs, forbidden_pairs, unpaired)
      prepare_energies(pairing_energy)
    end

    # @return [Integer] number of bases
    def length
      @sequence.length
    end

    # @return [Array<Array<Integer>>] forced pairs
    def forced_pairs
      @constraints.forced_pairs
    end

    # @return [Array<Array<Integer>>] forbidden pairs
    def forbidden_pairs
      @constraints.forbidden_pairs
    end

    # @return [Array<Integer>] bases forced to remain unpaired
    def unpaired_positions
      @constraints.unpaired_positions
    end

    # @return [Numeric, nil] caller supplied energy for pair (left, right)
    def pairing_energy(left, right)
      check_index(left.to_int)
      check_index(right.to_int)
      energy_at(left.to_int, right.to_int)
    end

    private

    def build_constraints(forced_pairs, forbidden_pairs, unpaired)
      Constraints.new(length, forced_pairs: forced_pairs,
                             forbidden_pairs: forbidden_pairs,
                             unpaired: unpaired)
    end
  end
end
