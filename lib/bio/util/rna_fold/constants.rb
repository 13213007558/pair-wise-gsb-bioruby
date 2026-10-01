# frozen_string_literal: true

# bio/util/rna_fold/constants.rb - Shared numerical constants
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Log weight of an impossible structure or substructure.
    NEG_INFINITY = -Float::INFINITY

    # Stable log(sum(exp(values))) in the shifted (max-subtracted) form.
    # @api private
    def self.log_sum_exp(terms)
      values = terms.reject { |value| value == NEG_INFINITY }
      return NEG_INFINITY if values.empty?
      maximum = values.max
      maximum + Math.log(values.sum { |value| Math.exp(value - maximum) })
    end
  end
end
