require 'bio/data/codontable'

module Bio
  module Alignment
    # Nei and Gojobori (1986) nonsynonymous/synonymous substitution distances.
    #
    # The calculator expects two aligned coding sequences with lengths that are
    # divisible by three. Gaps, ambiguous codons, stop codons, and replacement
    # paths that pass through stop codons can be excluded or reported as errors.
    class NeiGojobori
      MODES = [:exclude, :raise].freeze
      DEFAULT_GAP_REGEXP = /[-.]/.freeze
      UNAMBIGUOUS_CODON_REGEXP = /\A[acgt]{3}\z/.freeze
      BASES = %w[a c g t].freeze

      Result = Struct.new(:status, :status_reasons,
                          :nonsynonymous_sites, :synonymous_sites,
                          :nonsynonymous_nucleotide_differences,
                          :synonymous_nucleotide_differences,
                          :dn, :ds, :omega,
                          :excluded_columns,
                          :excluded_gap_columns,
                          :excluded_ambiguous_codon_columns,
                          :excluded_stop_codon_columns,
                          :comparable_codon_columns)

      attr_reader :codon_table, :gap_mode, :ambiguous_codon_mode,
                  :stop_codon_mode, :gap_regexp

      def initialize(codon_table: 1,
                     gap_mode: :exclude,
                     ambiguous_codon_mode: :exclude,
                     stop_codon_mode: :exclude,
                     gap_regexp: nil)
        @codon_table =
          if codon_table.is_a?(Bio::CodonTable)
            codon_table
          else
            Bio::CodonTable[codon_table]
          end
        @gap_mode = validate_mode(gap_mode, :gap_mode)
        @ambiguous_codon_mode = validate_mode(ambiguous_codon_mode, :ambiguous_codon_mode)
        @stop_codon_mode = validate_mode(stop_codon_mode, :stop_codon_mode)
        @gap_regexp = gap_regexp || DEFAULT_GAP_REGEXP
      end

      def calculate(seq1, seq2)
        codons1 = codons_from(seq1)
        codons2 = codons_from(seq2)
        if codons1.length != codons2.length
          raise ArgumentError, 'aligned coding sequences must have the same length'
        end

        n_sites = 0.0
        s_sites = 0.0
        n_differences = 0.0
        s_differences = 0.0
        comparable_columns = 0
        excluded = Hash.new(0)
        status_reasons = []

        codons1.each_with_index do |codon1, index|
          codon2 = codons2[index]
          reason = exclusion_reason(codon1, codon2)

          if reason
            excluded[reason] += 1
            excluded[:all] += 1
            handle_exclusion(reason, index)
            next
          end

          pair_n_sites, pair_s_sites = pair_codon_sites(codon1, codon2)
          nonsynonymous, synonymous, path_count =
            replacement_differences(codon1, codon2)

          if path_count.zero?
            excluded[:stop] += 1
            excluded[:all] += 1
            status_reasons << :no_valid_replacement_path
            handle_exclusion(:stop, index)
            next
          end

          n_sites += pair_n_sites
          s_sites += pair_s_sites
          n_differences += nonsynonymous
          s_differences += synonymous
          comparable_columns += 1
        end

        build_result(n_sites, s_sites, n_differences, s_differences,
                     comparable_columns, excluded, status_reasons)
      end

      def self.calculate(seq1, seq2, options = {})
        new(**options).calculate(seq1, seq2)
      end

      def self.dnds(seq1, seq2, options = {})
        calculate(seq1, seq2, options)
      end

      private

      def validate_mode(mode, option_name)
        return mode if MODES.include?(mode)

        raise ArgumentError,
              "#{option_name} must be one of #{MODES.map(&:inspect).join(', ')}"
      end

      def codons_from(sequence)
        string = sequence.to_s.downcase.tr('u', 't')
        if string.empty? || (string.length % 3).positive?
          raise ArgumentError, 'aligned coding sequences must have a length divisible by three'
        end

        string.scan(/.{3}/)
      end

      def exclusion_reason(codon1, codon2)
        [codon1, codon2].each do |codon|
          return :gap if gap_regexp.match(codon)
        end
        [codon1, codon2].each do |codon|
          return :ambiguous unless UNAMBIGUOUS_CODON_REGEXP.match(codon)
        end
        [codon1, codon2].each do |codon|
          return :stop if codon_table.stop_codon?(codon)
        end
        nil
      end

      def handle_exclusion(reason, index)
        mode = case reason
               when :gap then gap_mode
               when :ambiguous then ambiguous_codon_mode
               when :stop then stop_codon_mode
               end
        return unless mode == :raise

        labels = { gap: 'gap', ambiguous: 'ambiguous codon',
                   stop: 'stop codon' }
        label = labels.fetch(reason)
        raise ArgumentError, "#{label} found at aligned codon column #{index + 1}"
      end

      def pair_codon_sites(codon1, codon2)
        n1, s1 = codon_sites(codon1)
        n2, s2 = codon_sites(codon2)
        [(n1 + n2) / 2.0, (s1 + s2) / 2.0]
      end

      def codon_sites(codon)
        synonymous = 0
        amino_acid = codon_table[codon]

        3.times do |position|
          BASES.each do |base|
            next if codon[position] == base

            mutant = codon.dup
            mutant[position] = base
            mutant_amino_acid = codon_table[mutant]
            synonymous += 1 if mutant_amino_acid == amino_acid
          end
        end

        synonymous_sites = synonymous / 3.0
        [3.0 - synonymous_sites, synonymous_sites]
      end

      def replacement_differences(codon1, codon2)
        positions = differing_positions(codon1, codon2)
        return [0.0, 0.0, 1] if positions.empty?

        nonsynonymous = 0.0
        synonymous = 0.0
        valid_paths = 0

        positions.permutation do |path|
          current = codon1.dup
          previous_amino_acid = codon_table[codon1]
          path_nonsynonymous = 0
          path_synonymous = 0
          valid_path = true

          path.each do |position|
            current[position] = codon2[position]
            if codon_table.stop_codon?(current)
              valid_path = false
              break
            end

            amino_acid = codon_table[current]
            if amino_acid == previous_amino_acid
              path_synonymous += 1
            else
              path_nonsynonymous += 1
            end
            previous_amino_acid = amino_acid
          end

          next unless valid_path

          nonsynonymous += path_nonsynonymous
          synonymous += path_synonymous
          valid_paths += 1
        end

        return [0.0, 0.0, 0] if valid_paths.zero?

        [nonsynonymous / valid_paths, synonymous / valid_paths, valid_paths]
      end

      def differing_positions(codon1, codon2)
        (0...3).select { |position| codon1[position] != codon2[position] }
      end

      def build_result(n_sites, s_sites, n_differences, s_differences,
                       comparable_columns, excluded, status_reasons)
        pn = n_sites.positive? ? n_differences / n_sites : nil
        ps = s_sites.positive? ? s_differences / s_sites : nil
        dn, dn_reason = jukes_cantor_distance(pn, n_sites, true)
        ds, ds_reason = jukes_cantor_distance(ps, s_sites, false)

        status_reasons.concat([dn_reason, ds_reason].compact)
        status_reasons << :no_comparable_codon_columns if comparable_columns.zero?
        status_reasons.uniq!

        status = result_status(dn, ds, dn_reason, ds_reason, comparable_columns)
        omega = status == :ok ? dn / ds : nil

        Result.new(status, status_reasons,
                   n_sites, s_sites, n_differences, s_differences,
                   dn, ds, omega,
                   excluded[:all],
                   excluded[:gap],
                   excluded[:ambiguous],
                   excluded[:stop],
                   comparable_columns)
      end

      def jukes_cantor_distance(proportion, sites, nonsynonymous)
        if sites <= 0
          reason = if nonsynonymous
                     :no_nonsynonymous_sites
                   else
                     :no_synonymous_sites
                   end
          return [nil, reason]
        end
        return [0.0, nil] if proportion.nil? || proportion.zero?

        return [nil, :jukes_cantor_out_of_domain] if proportion >= 0.75

        [-0.75 * Math.log(1.0 - (4.0 / 3.0) * proportion), nil]
      end

      def result_status(dn, ds, dn_reason, ds_reason, comparable_columns)
        return :no_comparable_codon_columns if comparable_columns.zero?
        if dn_reason == :jukes_cantor_out_of_domain ||
           ds_reason == :jukes_cantor_out_of_domain
          return :jukes_cantor_out_of_domain
        end
        return :undefined_synonymous_distance if ds.nil?
        return :zero_synonymous_distance if ds.zero?
        return :undefined_nonsynonymous_distance if dn.nil?

        :ok
      end
    end
  end
end
