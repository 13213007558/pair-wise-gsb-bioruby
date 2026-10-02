# frozen_string_literal: true
#
# = bio/util/nei_gojobori.rb - Nei-Gojobori (1986) dN/dS estimation
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License::   The Ruby License
#

require 'bio/data/codontable'

module Bio

  # = Bio::NeiGojobori
  #
  # Estimates the numbers of synonymous and nonsynonymous substitutions
  # per site (dS and dN) and their ratio (omega) between two aligned
  # coding sequences, using the method of Nei and Gojobori (1986) with
  # the Jukes-Cantor correction for multiple hits.
  #
  # The calculator works directly on a pair of aligned nucleotide
  # sequences (String, Bio::Sequence::NA, or any object whose +to_s+
  # gives the aligned sequence) or on a Bio::Alignment object holding
  # exactly two sequences.  Translation and stop-codon detection are
  # delegated to Bio::CodonTable, so any genetic code supported by
  # Bio::CodonTable can be used.
  #
  # Columns containing gaps, ambiguous (degenerate) nucleotides, or
  # stop codons are handled according to explicit per-category
  # policies (:exclude to skip the column, :raise to raise
  # Bio::NeiGojobori::Error).  Codon pairs differing at more than one
  # position are resolved by averaging over all mutational paths that
  # do not pass through a stop codon, instead of fixing a single
  # arbitrary order of changes.
  #
  # == Usage
  #
  #   require 'bio'
  #
  #   result = Bio::NeiGojobori.run(seq1, seq2,
  #                                 codon_table: Bio::CodonTable[1],
  #                                 gap: :exclude,
  #                                 ambiguous: :exclude,
  #                                 stop: :exclude)
  #   result.status                  # => :ok, :zero_ds, :saturated, or :no_sites
  #   result.dN                      # => Float or nil
  #   result.dS                      # => Float or nil
  #   result.omega                   # => Float or nil (nil unless status == :ok)
  #   result.synonymous_sites        # => Float
  #   result.nonsynonymous_sites     # => Float
  #   result.synonymous_differences  # => Float
  #   result.nonsynonymous_differences # => Float
  #   result.compared_columns        # => Integer
  #   result.excluded_columns        # => Integer
  #
  # From a Bio::Alignment object with exactly two sequences:
  #
  #   result = Bio::NeiGojobori.run_alignment(alignment)
  #
  # When the Jukes-Cantor correction falls outside its domain
  # (proportion of differences >= 3/4) the corresponding distance is
  # reported as nil and the status becomes :saturated.  When dS is
  # zero the ratio is undefined: omega is nil and the status becomes
  # :zero_ds.  In neither case is a normal-looking ratio returned.
  #
  # == References
  #
  # * Nei M. and Gojobori T.  Simple methods for estimating the numbers
  #   of synonymous and nonsynonymous nucleotide substitutions.
  #   Mol Biol Evol. 1986 3(5): 418-426.
  #
  class NeiGojobori

    # Error raised when the input sequences cannot be analyzed, or when
    # a handling policy of :raise is triggered.
    class Error < RuntimeError; end

    # Result of a Nei-Gojobori calculation.  Distances and the dN/dS
    # ratio are nil when they cannot be defined; always check #status
    # (or #ok?) before interpreting them.
    class Result
      # Computation status: :ok, :zero_ds, :saturated, or :no_sites.
      attr_reader :status
      # Number of synonymous sites (averaged over the two sequences).
      attr_reader :synonymous_sites
      # Number of nonsynonymous sites (averaged over the two sequences).
      attr_reader :nonsynonymous_sites
      # Number of synonymous differences (averaged over valid paths).
      attr_reader :synonymous_differences
      # Number of nonsynonymous differences (averaged over valid paths).
      attr_reader :nonsynonymous_differences
      # Proportion of synonymous differences (synonymous_differences /
      # synonymous_sites), or nil when synonymous_sites is zero.
      attr_reader :pS
      # Proportion of nonsynonymous differences, or nil.
      attr_reader :pN
      # Jukes-Cantor corrected synonymous distance, or nil.
      attr_reader :dS
      # Jukes-Cantor corrected nonsynonymous distance, or nil.
      attr_reader :dN
      # dN/dS ratio, or nil unless status is :ok.
      attr_reader :omega
      # Number of codon columns used in the calculation.
      attr_reader :compared_columns
      # Number of codon columns excluded from the calculation.
      attr_reader :excluded_columns

      def initialize(values) #:nodoc:
        values.each do |key, value|
          instance_variable_set("@#{key}", value)
        end
      end

      # True when dN, dS and omega are all defined.
      def ok?
        @status == :ok
      end

      # True when the Jukes-Cantor correction fell outside its domain
      # for dN and/or dS.
      def saturated?
        @status == :saturated
      end

      # dN/dS ratio (alias of omega).
      alias dN_dS omega
    end

    # Handling policies for gap, ambiguous, and stop-codon columns.
    HANDLING_POLICIES = %i[exclude raise].freeze

    # Nucleotides used when counting synonymous sites and paths.
    BASES = %w[a c g t].freeze

    # Upper limit of the Jukes-Cantor correction domain.
    JC_SATURATION_LIMIT = 0.75

    # Creates a new calculator.
    #
    # Options:
    # [codon_table] Bio::CodonTable used for translation and stop-codon
    #               detection (default: Bio::CodonTable[1]).
    # [gap]         policy for columns containing gap characters in
    #               either sequence: :exclude (default) or :raise.
    # [ambiguous]   policy for columns containing degenerate nucleotide
    #               codes (n, r, y, ...): :exclude (default) or :raise.
    # [stop]        policy for columns where either codon is a stop
    #               codon: :exclude (default) or :raise.
    def initialize(codon_table: nil, gap: :exclude, ambiguous: :exclude,
                   stop: :exclude)
      @codon_table = codon_table || Bio::CodonTable[1]
      @gap_policy = check_policy(:gap, gap)
      @ambiguous_policy = check_policy(:ambiguous, ambiguous)
      @stop_policy = check_policy(:stop, stop)
    end

    # Runs the calculation on a pair of aligned sequences.
    #
    #   Bio::NeiGojobori.run(seq1, seq2, stop: :raise)
    def self.run(seq1, seq2, **options)
      new(**options).run(seq1, seq2)
    end

    # Runs the calculation on a Bio::Alignment (or any object
    # responding to each_seq) holding exactly two sequences.
    #
    #   Bio::NeiGojobori.run_alignment(alignment)
    def self.run_alignment(alignment, **options)
      new(**options).run_alignment(alignment)
    end

    # Runs the calculation on a pair of aligned sequences and returns
    # a Bio::NeiGojobori::Result.
    def run(seq1, seq2)
      str1 = normalize_sequence(seq1)
      str2 = normalize_sequence(seq2)
      unless str1.length == str2.length
        raise Error, "aligned sequences must have the same length (#{str1.length} != #{str2.length})"
      end
      unless (str1.length % 3).zero?
        raise Error, "aligned sequence length #{str1.length} is not a multiple of 3"
      end

      s_sites = 0.0
      n_sites = 0.0
      s_diffs = 0.0
      n_diffs = 0.0
      compared = 0
      excluded = 0

      (str1.length / 3).times do |k|
        codon1 = str1[k * 3, 3]
        codon2 = str2[k * 3, 3]
        if exclude_column?(codon1, codon2)
          excluded += 1
          next
        end
        diffs = codon_differences(codon1, codon2)
        unless diffs
          # No mutational path avoids a stop codon.
          excluded += 1
          next
        end
        compared += 1
        syn1 = synonymous_sites(codon1)
        syn2 = synonymous_sites(codon2)
        s_sites += (syn1 + syn2) / 2.0
        n_sites += (6.0 - syn1 - syn2) / 2.0
        s_diffs += diffs[0]
        n_diffs += diffs[1]
      end

      build_result(s_sites, n_sites, s_diffs, n_diffs, compared, excluded)
    end

    # Runs the calculation on an alignment holding exactly two
    # sequences and returns a Bio::NeiGojobori::Result.
    def run_alignment(alignment)
      seqs = []
      alignment.each_seq { |seq| seqs << seq }
      unless seqs.size == 2
        raise Error, "alignment must contain exactly 2 sequences (got #{seqs.size})"
      end
      run(seqs[0], seqs[1])
    end

    private

    def check_policy(name, policy)
      return policy if HANDLING_POLICIES.include?(policy)

      policies = HANDLING_POLICIES.map(&:inspect).join(', ')
      raise ArgumentError,
            "invalid #{name} policy #{policy.inspect} (expected one of #{policies})"
    end

    def normalize_sequence(seq)
      seq.to_s.downcase.tr('u', 't')
    end

    # Returns true when the column must be skipped.  Raises
    # Bio::NeiGojobori::Error when the relevant policy is :raise.
    def exclude_column?(codon1, codon2)
      [codon1, codon2].each do |codon|
        if codon =~ /[^a-z]/
          handle_policy(@gap_policy, :gap,
                        "gap character in codon '#{codon}'")
          return true
        elsif codon =~ /[^acgt]/
          handle_policy(@ambiguous_policy, :ambiguous,
                        "ambiguous nucleotide in codon '#{codon}'")
          return true
        elsif @codon_table.stop_codon?(codon)
          handle_policy(@stop_policy, :stop,
                        "stop codon '#{codon}'")
          return true
        end
      end
      false
    end

    def handle_policy(policy, category, message)
      return unless policy == :raise

      raise Error, "#{category} column found: #{message}"
    end

    # Number of synonymous sites of a codon (Nei-Gojobori definition):
    # for each position, one third of the number of synonymous changes
    # among the three possible nucleotide substitutions.
    def synonymous_sites(codon)
      aa = @codon_table[codon]
      sites = 0.0
      3.times do |pos|
        synonymous = BASES.count do |base|
          next false if base == codon[pos]

          @codon_table[replace_base(codon, pos, base)] == aa
        end
        sites += synonymous / 3.0
      end
      sites
    end

    # Synonymous and nonsynonymous differences between two codons.
    # When the codons differ at more than one position, all mutational
    # orders are considered and the counts are averaged over the paths
    # that do not pass through a stop codon.  Returns nil when no such
    # path exists.
    def codon_differences(codon1, codon2)
      differing = (0..2).select { |pos| codon1[pos] != codon2[pos] }
      case differing.size
      when 0
        [0.0, 0.0]
      when 1
        if @codon_table[codon1] == @codon_table[codon2]
          [1.0, 0.0]
        else
          [0.0, 1.0]
        end
      else
        averaged_path_differences(codon1, codon2, differing)
      end
    end

    def averaged_path_differences(codon1, codon2, differing)
      paths = differing.permutation.filter_map do |order|
        path_difference_counts(codon1, codon2, order)
      end
      return nil if paths.empty?

      s_total = paths.sum { |counts| counts[0] }
      n_total = paths.sum { |counts| counts[1] }
      [s_total.fdiv(paths.size), n_total.fdiv(paths.size)]
    end

    # Synonymous and nonsynonymous change counts along one mutational
    # order, or nil when the path passes through a stop codon.
    def path_difference_counts(codon1, codon2, order)
      synonymous = 0.0
      current = codon1
      order.each do |pos|
        nxt = replace_base(current, pos, codon2[pos])
        return nil if @codon_table.stop_codon?(nxt)

        synonymous += 1 if @codon_table[current] == @codon_table[nxt]
        current = nxt
      end
      [synonymous, order.size - synonymous]
    end

    def replace_base(codon, pos, base)
      replaced = codon.dup
      replaced[pos] = base
      replaced
    end

    # Jukes-Cantor correction; nil when p is outside the domain.
    def jukes_cantor(p)
      return nil if p >= JC_SATURATION_LIMIT
      return 0.0 if p.zero?

      -0.75 * Math.log(1.0 - (4.0 * p / 3.0))
    end

    def build_result(s_sites, n_sites, s_diffs, n_diffs, compared, excluded)
      ps = s_sites.positive? ? s_diffs / s_sites : nil
      pn = n_sites.positive? ? n_diffs / n_sites : nil
      ds = ps && jukes_cantor(ps)
      dn = pn && jukes_cantor(pn)
      status, omega = resolve_status(ps, pn, ds, dn)
      Result.new(status: status,
                 synonymous_sites: s_sites,
                 nonsynonymous_sites: n_sites,
                 synonymous_differences: s_diffs,
                 nonsynonymous_differences: n_diffs,
                 pS: ps, pN: pn, dS: ds, dN: dn, omega: omega,
                 compared_columns: compared,
                 excluded_columns: excluded)
    end

    def resolve_status(ps, pn, ds, dn)
      return [:no_sites, nil] if ps.nil? || pn.nil?
      return [:saturated, nil] if ds.nil? || dn.nil?
      return [:zero_ds, nil] if ds.zero?

      [:ok, dn / ds]
    end
  end
end
