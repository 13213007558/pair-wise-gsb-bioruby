# frozen_string_literal: true
#
# = bio/util/pcr.rb - Sequence-level in-silico PCR simulation
#
# == Bio::PCR - Simulate PCR amplification against a template sequence.
#
# Given a pair of primers (which may contain IUPAC degenerate bases),
# this class searches a template sequence for primer binding sites on
# both strands and reports all possible amplicons.
#
# == Example
#
#  template = Bio::Sequence::NA.new('ACGT' * 250)
#  pcr = Bio::PCR.new(template)
#  products = pcr.amplify('ACGTACGTAC', 'GTACGTACGT',
#                         max_mismatches: 1,
#                         max_three_prime_mismatches: 0,
#                         min_product_length: 50,
#                         max_product_length: 500)
#  products.each do |product|
#    puts [product.start, product.stop, product.length].join("\t")
#    puts product.sequence
#    p product.forward_match.mismatches
#    p product.reverse_match.mismatches
#  end
#

require 'bio/sequence'

module Bio
  # = Bio::PCR
  # Sequence-level in-silico PCR simulation with support for IUPAC
  # degenerate bases in primers, mismatch tolerance, stricter 3'-end
  # mismatch rules, product length limits and circular templates.
  class PCR
    # IUPAC degenerate base symbols mapped to the bases they represent.
    IUPAC_BASES = {
      'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T', 'U' => 'T',
      'R' => 'AG', 'Y' => 'CT', 'S' => 'GC', 'W' => 'AT',
      'K' => 'GT', 'M' => 'AC',
      'B' => 'CGT', 'D' => 'AGT', 'H' => 'ACT', 'V' => 'ACG',
      'N' => 'ACGT'
    }.freeze

    # A single primer binding site found on the template.
    class PrimerMatch
      # 0-based start position of the site on the plus strand of the
      # template. For circular templates the position is always
      # smaller than the template length.
      attr_reader :start

      # 0-based exclusive end position of the site on the plus strand.
      # For circular templates the value may exceed the template
      # length when the site wraps around the origin.
      attr_reader :stop

      # Strand the primer binds to, as a Symbol: +:plus+ when the
      # primer sequence is identical to the plus strand (right-facing
      # site), +:minus+ when it is identical to the minus strand
      # (left-facing site).
      attr_reader :strand

      # Mismatch positions as an Array of 0-based offsets counted
      # from the 5' end of the primer.
      attr_reader :mismatches

      def initialize(start, stop, strand, mismatches) # :nodoc:
        @start = start
        @stop = stop
        @strand = strand
        @mismatches = mismatches.freeze
      end

      # Length of the binding site (i.e. the primer length).
      def length
        @stop - @start
      end

      # Number of mismatches located within the last +window+ bases
      # at the 3' end of the primer.
      def three_prime_mismatches(window)
        border = length - window
        @mismatches.count { |pos| pos >= border }
      end
    end

    # A predicted PCR product (amplicon).
    class Product
      # 0-based start position of the amplicon on the plus strand of
      # the template (smaller than the template length).
      attr_reader :start

      # 0-based exclusive end position of the amplicon on the plus
      # strand. For circular templates the value may exceed the
      # template length when the amplicon wraps around the origin
      # (at most one full circle).
      attr_reader :stop

      # Amplicon sequence as a Bio::Sequence::NA, oriented from the
      # forward primer to the reverse primer.
      attr_reader :sequence

      # PrimerMatch of the forward primer.
      attr_reader :forward_match

      # PrimerMatch of the reverse primer.
      attr_reader :reverse_match

      def initialize(start, stop, sequence, forward_match, reverse_match) # :nodoc:
        @start = start
        @stop = stop
        @sequence = sequence
        @forward_match = forward_match
        @reverse_match = reverse_match
      end

      # Orientation of the amplicon relative to the template plus
      # strand, as a Symbol: +:plus+ or +:minus+. The amplicon
      # sequence is always reported 5'->3' from the forward primer
      # to the reverse primer.
      def strand
        @forward_match.strand
      end

      # Length of the amplicon in base pairs.
      def length
        @stop - @start
      end
    end

    # Matching rules and product length limits of a simulation.
    class Rules
      # Maximum total number of mismatches allowed per primer.
      attr_reader :max_mismatches

      # Number of bases at the 3' end of each primer checked with the
      # stricter rule.
      attr_reader :three_prime_window

      # Maximum number of mismatches allowed inside the 3'-end window.
      attr_reader :max_three_prime_mismatches

      # Minimum amplicon length in bp.
      attr_reader :min_product_length

      # Maximum amplicon length in bp.
      attr_reader :max_product_length

      def initialize(max_mismatches: 0, three_prime_window: 5, # :nodoc:
                     max_three_prime_mismatches: 0,
                     min_product_length: 0, max_product_length: nil)
        @max_mismatches = max_mismatches
        @three_prime_window = three_prime_window
        @max_three_prime_mismatches = max_three_prime_mismatches
        @min_product_length = min_product_length
        @max_product_length = max_product_length
      end
    end

    # Option keys accepted by Bio::PCR#amplify.
    RULE_KEYS = %i[max_mismatches three_prime_window
                   max_three_prime_mismatches
                   min_product_length max_product_length].freeze

    # Template sequence as a Bio::Sequence::NA.
    attr_reader :template

    # Whether the template is treated as a circular molecule.
    attr_reader :circular

    # Creates a new simulator. +template+ may be a String or a
    # Bio::Sequence::NA. When +circular+ is true, binding sites and
    # products may wrap around the origin, but products are limited
    # to at most one full circle of the template.
    def initialize(template, circular: false)
      seq = template.to_s.upcase
      raise ArgumentError, 'template sequence is empty' if seq.empty?

      @sequence = seq
      @template = Bio::Sequence::NA.new(seq)
      @circular = circular
    end

    # Simulates PCR with the given primer pair and returns an Array
    # of Bio::PCR::Product, sorted by position. Every distinct pair
    # of binding sites yields its own product.
    #
    # Options:
    # [+max_mismatches+]
    #   Maximum total number of mismatches allowed per primer
    #   (default: 0, perfect matches only).
    # [+three_prime_window+]
    #   Number of bases at the 3' end of each primer that follow the
    #   stricter rule below (default: 5).
    # [+max_three_prime_mismatches+]
    #   Maximum number of mismatches allowed inside the 3'-end window
    #   (default: 0).
    # [+min_product_length+]
    #   Minimum amplicon length in bp (default: 0).
    # [+max_product_length+]
    #   Maximum amplicon length in bp (default: the template length;
    #   circular templates never exceed one full circle).
    def amplify(forward_primer, reverse_primer, **options)
      rules = build_rules(options)
      fwd = normalize_primer(forward_primer)
      rev = normalize_primer(reverse_primer)
      scan = @circular ? @sequence * 2 : @sequence
      products = collect_products(scan, fwd, rev, rules)
      products.uniq! { |product| [product.start, product.stop, product.strand] }
      products.sort_by! { |product| [product.start, product.stop] }
      products
    end

    private

    def build_rules(options)
      unknown = options.keys - RULE_KEYS
      raise ArgumentError, "unknown options: #{unknown.join(', ')}" unless unknown.empty?

      Rules.new(**options)
    end

    def normalize_primer(primer)
      seq = primer.to_s.upcase
      raise ArgumentError, 'primer sequence is empty' if seq.empty?

      bad = seq.each_char.find { |char| !IUPAC_BASES.key?(char) }
      raise ArgumentError, "invalid base in primer: #{bad}" if bad

      seq
    end

    def effective_max_length(rules)
      limit = @template.length
      max = rules.max_product_length
      return limit if max.nil? || max > limit

      max
    end

    def collect_products(scan, fwd, rev, rules)
      products = []
      [[fwd, rev, :plus], [rev, fwd, :minus]].each do |right, left, strand|
        products.concat(pair_sites(scan, right, left, strand, rules))
      end
      products
    end

    def pair_sites(scan, right, left, strand, rules)
      b_sites = left_sites(scan, left, rules)
      b_sites += shifted_sites(b_sites) if @circular
      products = []
      right_sites(scan, right, rules).each do |site_a|
        b_sites.each do |site_b|
          product = build_product(scan, site_a, site_b, strand, rules)
          products << product if product
        end
      end
      products
    end

    # Copies of the given sites shifted by one template length, so
    # that products wrapping around the origin of a circular template
    # (at most one full circle) can be formed.
    def shifted_sites(sites)
      length = @template.length
      sites.map do |site|
        PrimerMatch.new(site.start + length, site.stop + length,
                        site.strand, site.mismatches)
      end
    end

    def build_product(scan, site_a, site_b, strand, rules)
      return nil unless site_a.start < site_b.start

      length = site_b.stop - site_a.start
      return nil if length < rules.min_product_length
      return nil if length > effective_max_length(rules)

      region = scan[site_a.start...site_b.stop]
      if strand == :plus
        Product.new(site_a.start, site_b.stop,
                    Bio::Sequence::NA.new(region), site_a, site_b)
      else
        Product.new(site_a.start, site_b.stop,
                    Bio::Sequence::NA.new(reverse_complement(region)),
                    site_b, site_a)
      end
    end

    # Sites where the primer itself matches the plus strand
    # (right-facing sites).
    def right_sites(scan, primer, rules)
      scan_sites(scan, primer, :end, rules).map do |start, mismatches|
        PrimerMatch.new(start, start + primer.length, :plus, mismatches)
      end
    end

    # Sites where the reverse complement of the primer matches the
    # plus strand (left-facing sites). Mismatch offsets are converted
    # back to 5'->3' offsets of the primer itself. Note that the 3'
    # end of the primer corresponds to the 5' end of its reverse
    # complement, so the strict 3'-end window is located at the
    # beginning of the scanned pattern.
    def left_sites(scan, primer, rules)
      rc = reverse_complement(primer)
      scan_sites(scan, rc, :start, rules).map do |start, mismatches|
        converted = mismatches.map { |pos| primer.length - 1 - pos }
        PrimerMatch.new(start, start + primer.length, :minus, converted)
      end
    end

    # Returns an Array of [start, mismatch_offsets] for all binding
    # sites of +pattern+ on +scan+ satisfying the mismatch rules.
    # +three_prime_at+ is +:end+ when the 3' end of the primer is at
    # the end of +pattern+, or +:start+ when it is at the beginning
    # (i.e. +pattern+ is the reverse complement of the primer).
    def scan_sites(scan, pattern, three_prime_at, rules)
      sites = []
      each_site_start(scan, pattern.length) do |start|
        mismatches = mismatches_at(scan, start, pattern)
        next if mismatches.nil?
        next if mismatches.size > rules.max_mismatches
        next if three_prime_exceeded?(mismatches, pattern.length,
                                      three_prime_at, rules)

        sites << [start, mismatches]
      end
      sites
    end

    def three_prime_exceeded?(mismatches, pattern_length, three_prime_at, rules)
      window = rules.three_prime_window
      count = if three_prime_at == :start
                mismatches.count { |pos| pos < window }
              else
                mismatches.count { |pos| pos >= pattern_length - window }
              end
      count > rules.max_three_prime_mismatches
    end

    def each_site_start(scan, pattern_length)
      if @circular
        return if pattern_length > @template.length

        (0...@template.length).each { |start| yield start }
      else
        max_start = scan.length - pattern_length
        (0..max_start).each { |start| yield start } unless max_start.negative?
      end
    end

    # Returns an Array of mismatch offsets (0-based from the 5' end of
    # the pattern), or nil when the site extends beyond the scan
    # sequence.
    def mismatches_at(scan, start, pattern)
      mismatches = []
      pattern.length.times do |offset|
        base = scan[start + offset]
        return nil if base.nil?

        mismatches << offset unless base_matches?(pattern[offset], base)
      end
      mismatches
    end

    def base_matches?(pattern_base, template_base)
      pattern_set = IUPAC_BASES[pattern_base]
      template_set = IUPAC_BASES[template_base]
      return false if pattern_set.nil? || template_set.nil?

      pattern_set.each_char.any? { |base| template_set.include?(base) }
    end

    def reverse_complement(sequence)
      sequence.tr('ACGTU', 'TGCAA').reverse
    end
  end
end
