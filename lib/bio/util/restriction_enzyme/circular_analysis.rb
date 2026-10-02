#
# bio/util/restriction_enzyme/circular_analysis.rb - Circular and methylation-aware restriction analysis
#
# Author::    BioRuby contributors
# Copyright:: Copyright (c) 2026 BioRuby project
# License::   The Ruby License
#

require 'bio/sequence'

module Bio

require 'bio/util/restriction_enzyme' unless const_defined?(:RestrictionEnzyme)

class RestrictionEnzyme

  # = Description
  #
  # Bio::RestrictionEnzyme::CircularAnalysis performs restriction
  # digestion analysis on linear or circular (e.g. plasmid) DNA.
  #
  # Unlike Bio::RestrictionEnzyme::Analysis it does not rely on the
  # built-in REBASE enzyme database.  Instead, the caller supplies
  # explicit enzyme rules: the recognition site, where each strand is
  # cut, and how methylation affects the enzyme.  Recognition sites
  # may cross the first base of a circular sequence (the origin).
  #
  # The result reports:
  #
  # * +cuts+          - cut (nick) positions on each strand
  # * +nicks+         - single-strand cuts
  # * +fragments+     - the resulting fragments with their sticky ends
  # * +blocked_sites+ - recognition sites that were shielded by methylation
  #
  # A circular molecule with exactly one double-strand cut yields a
  # single linear fragment.  Rotating a circular input sequence yields
  # an equivalent result, shifted only in coordinates.
  #
  # = Usage
  #
  #   require 'bio'
  #   require 'bio/util/restriction_enzyme/circular_analysis'
  #
  #   eco_ri = { name: 'EcoRI', site: 'g^aattc', cut_complement: 5,
  #              methylation: :blocked }
  #
  #   result = Bio::RestrictionEnzyme::CircularAnalysis.cut(
  #     'atgcgaattcgc',
  #     eco_ri,
  #     circular: true,
  #     methylation: { primary: [5], complement: [] }
  #   )
  #
  #   result.cuts          # => { primary: [...], complement: [...] }
  #   result.fragments     # => array of Fragment
  #   result.blocked_sites # => array of Site
  #
  # = Coordinates
  #
  # All positions are 0-based offsets into the primary (given) strand.
  # A cut position +k+ means the strand is broken immediately before
  # the base at index +k+.  Complement-strand positions use the same
  # coordinate system: complement position +i+ pairs with primary
  # position +i+.  For circular sequences all coordinates are taken
  # modulo the sequence length.
  #
  class CircularAnalysis

    # A caller-supplied enzyme rule.
    #
    # * +name+ - enzyme name (defaults to the site pattern)
    # * +site+ - recognition site on the primary strand.  May contain
    #   IUPAC ambiguity codes.  A caret (+^+) marks the primary-strand
    #   cut position.
    # * +cut_primary+ - offset from the site start to the primary
    #   strand cut (the strand is cut before the base at this offset).
    #   Overrides a caret in +site+.  Pass +false+ for enzymes that do
    #   not cut the primary strand.
    # * +cut_complement+ - offset to the complement-strand cut, in
    #   primary-strand coordinates.  Defaults to +cut_primary+ (blunt
    #   cut).  Pass +false+ for nicking enzymes that leave the
    #   complement strand intact.
    # * +methylation+ - sensitivity rule:
    #   * +:none+ (default) - methylation never blocks cutting
    #   * +:blocked+ - any methylated base inside the recognition
    #     site, on either strand, blocks cutting
    #   * +:blocked_primary+ / +:blocked_complement+ - only
    #     methylation of the given strand blocks cutting
    #   * any callable - called with a Site object; a truthy return
    #     value means cutting is blocked
    class EnzymeRule
      attr_reader :name, :site, :cut_primary, :cut_complement, :methylation

      def initialize(name: nil, site: nil, cut_primary: nil,
                     cut_complement: nil, methylation: :none)
        raise ArgumentError, 'site is required' if site.nil? || site.to_s.empty?
        site = site.to_s.downcase
        if site.include?('^')
          cut_primary = site.index('^') if cut_primary.nil?
          site = site.delete('^')
        end
        raise ArgumentError, 'cut position is required (use a caret in site or cut_primary:)' if cut_primary.nil?
        @name = (name || site).to_s
        @site = site
        @cut_primary = cut_primary
        @cut_complement = cut_complement.nil? ? cut_primary : cut_complement
        @methylation = methylation || :none
      end

      # Length of the recognition site in bases.
      def site_length
        @site.length
      end
    end # EnzymeRule

    # A recognition site found on the sequence.  Also used to report
    # sites whose cutting was blocked by methylation, and passed to
    # custom methylation rule callables.
    #
    # For circular sequences a site may cross the origin; in that case
    # +stop+ is greater than the sequence length.
    Site = Struct.new(:enzyme, :start, :stop, :sequence,
                      :primary_methylation, :complement_methylation)

    # A single-strand cut.  +strand+ is +:primary+ or +:complement+,
    # +position+ is the index of the base immediately after the cut,
    # and +enzyme+ is the name of the enzyme rule responsible.
    Nick = Struct.new(:strand, :position, :enzyme)

    # A fragment end.  +type+ is +:five_prime+ (5' overhang),
    # +:three_prime+ (3' overhang) or +:blunt+.  +overhang+ holds the
    # overhang sequence (empty for blunt ends).
    End = Struct.new(:type, :overhang)

    # A double-stranded fragment.  +primary+ and +complement+ are
    # padded with spaces so the strands align as in the original
    # molecule.  The +complement+ strand is written left to right in
    # its 3' -> 5' direction, base-paired with +primary+.
    # +left_end+ and +right_end+ are End objects.
    Fragment = Struct.new(:primary, :complement, :left_end, :right_end)

    # The result of a digestion.
    class Result
      # Array of Fragment objects.
      attr_reader :fragments
      # Hash with +:primary+ and +:complement+ keys, each holding a
      # sorted array of unique cut positions.
      attr_reader :cuts
      # Array of Nick objects (single-strand cuts).
      attr_reader :nicks
      # Array of Site objects that matched an enzyme rule but were
      # shielded from cutting by methylation.
      attr_reader :blocked_sites

      def initialize(fragments, cuts, nicks, blocked_sites)
        @fragments = fragments
        @cuts = cuts
        @nicks = nicks
        @blocked_sites = blocked_sites
      end
    end # Result

    # Convenience class method.  See #cut.
    def self.cut(sequence, *rules, circular: true, methylation: nil)
      new(sequence, circular: circular, methylation: methylation).cut(*rules)
    end

    # Normalizes a user-supplied rule into an EnzymeRule.  Accepts an
    # EnzymeRule, a Hash of EnzymeRule keywords, or a String site
    # pattern such as 'g^aattc' (blunt complement cut at the caret).
    def self.normalize_rule(rule)
      case rule
      when EnzymeRule
        rule
      when Hash
        EnzymeRule.new(**rule)
      when String
        EnzymeRule.new(site: rule)
      else
        raise ArgumentError, "cannot use #{rule.inspect} as an enzyme rule"
      end
    end

    # ---
    # *Arguments*
    # * +sequence+: the DNA sequence (primary strand, 5' -> 3')
    # * +circular+: +true+ if the molecule is circular (default)
    # * +methylation+: +Hash+ with +:primary+ and +:complement+ keys,
    #   each an array of methylated base positions (0-based indices
    #   into the primary strand)
    def initialize(sequence, circular: true, methylation: nil)
      @sequence = sequence.to_s.downcase
      @circular = circular ? true : false
      methylation ||= {}
      @methylation = {
        primary: methylation_positions(methylation, :primary),
        complement: methylation_positions(methylation, :complement)
      }
    end

    # Digests the sequence with the given enzyme rules.
    #
    # ---
    # *Arguments*
    # * +rules+: enzyme rules (EnzymeRule, Hash or String, see
    #   CircularAnalysis.normalize_rule)
    # *Returns*:: CircularAnalysis::Result
    def cut(*rules)
      rules = rules.flatten.map { |rule| self.class.normalize_rule(rule) }

      if @sequence.empty?
        return Result.new([], { primary: [], complement: [] }, [], [])
      end

      actions = [] # [rule, start] pairs that cut
      blocked = []
      rules.each do |rule|
        find_match_locations(rule).each do |start|
          site = build_site(rule, start)
          if blocked_by_methylation?(rule, site)
            blocked << site
          else
            actions << [rule, start]
          end
        end
      end

      breaks = [] # [primary_cut, complement_cut] pairs
      nicks = []
      primary_cut_positions = []
      complement_cut_positions = []
      actions.each do |rule, start|
        p_cut = cut_position(start, rule.cut_primary)
        c_cut = cut_position(start, rule.cut_complement)
        primary_cut_positions << p_cut if p_cut
        complement_cut_positions << c_cut if c_cut
        if p_cut && c_cut
          breaks << [p_cut, c_cut]
        elsif p_cut
          nicks << Nick.new(:primary, p_cut, rule.name)
        elsif c_cut
          nicks << Nick.new(:complement, c_cut, rule.name)
        end
      end

      cuts = {
        primary: primary_cut_positions.uniq.sort,
        complement: complement_cut_positions.uniq.sort
      }
      Result.new(build_fragments(breaks.uniq.sort), cuts, nicks, blocked)
    end

    #######
    private
    #######

    # Normalizes the methylated positions for one strand.
    def methylation_positions(methylation, strand)
      return [] if @sequence.empty?
      positions = methylation[strand] || methylation[strand.to_s] || []
      positions.map do |pos|
        pos = Integer(pos)
        if @circular
          pos % seq_length
        elsif pos < 0 || pos >= seq_length
          raise ArgumentError, "methylation position #{pos} is out of range"
        else
          pos
        end
      end.uniq.sort
    end

    # Returns the sorted match start positions of the rule's
    # recognition site.  For circular sequences the site may cross
    # the origin.
    def find_match_locations(rule)
      len = rule.site_length
      text = @circular ? @sequence * (len / seq_length.to_f).ceil + @sequence : @sequence
      re = Bio::Sequence::NA.new(rule.site).to_re
      locations = []
      index = 0
      while (md = text.match(re, index)) && md.begin(0) < seq_length
        locations << md.begin(0)
        index = md.begin(0) + 1
      end
      locations
    end

    # Builds a Site describing the recognition site starting at
    # +start+, including the methylated bases inside it.
    def build_site(rule, start)
      len = rule.site_length
      positions = (start...(start + len)).map { |i| i % seq_length }
      Site.new(
        rule.name,
        start,
        start + len,
        @circular ? circular_subseq(start, len) : @sequence[start, len],
        positions & @methylation[:primary],
        positions & @methylation[:complement]
      )
    end

    # Returns true if the rule's methylation sensitivity blocks
    # cutting at the given site.
    def blocked_by_methylation?(rule, site)
      case rule.methylation
      when :none, :insensitive
        false
      when :blocked
        !(site.primary_methylation.empty? && site.complement_methylation.empty?)
      when :blocked_primary
        !site.primary_methylation.empty?
      when :blocked_complement
        !site.complement_methylation.empty?
      else
        if rule.methylation.respond_to?(:call)
          rule.methylation.call(site) ? true : false
        else
          raise ArgumentError, "unknown methylation rule: #{rule.methylation.inspect}"
        end
      end
    end

    # Translates a cut offset within a recognition site into an
    # absolute cut position.  Returns nil when the strand is not cut.
    def cut_position(start, offset)
      return nil if offset.nil? || offset == false
      position = start + offset
      @circular ? position % seq_length : position
    end

    # Builds the fragments from the double-strand breaks.
    # +breaks+ is a sorted array of unique [primary_cut,
    # complement_cut] pairs.
    def build_fragments(breaks)
      if @circular
        return [] if breaks.empty?
        breaks.each_with_index.map do |current, i|
          nxt = breaks[(i + 1) % breaks.size]
          build_fragment(current, nxt, breaks.size == 1)
        end
      else
        bounds = [[0, 0]] + breaks + [[seq_length, seq_length]]
        bounds.each_cons(2).map do |left, right|
          build_fragment(left, right, false)
        end.reject { |fragment| fragment.primary.strip.empty? && fragment.complement.strip.empty? }
      end
    end

    # Builds one fragment spanning from break +left+ to break +right+.
    def build_fragment(left, right, full_circle)
      p_left, c_left = left
      p_right, c_right = right
      primary = strand_subseq(p_left, p_right, full_circle)
      complement = forward_complement(strand_subseq(c_left, c_right, full_circle))

      left_shift = signed_diff(c_left, p_left)
      right_shift = signed_diff(c_right, p_right)
      if left_shift > 0
        complement = (' ' * left_shift) + complement
      elsif left_shift < 0
        primary = (' ' * (-left_shift)) + primary
      end
      if right_shift > 0
        primary += ' ' * right_shift
      elsif right_shift < 0
        complement += ' ' * (-right_shift)
      end

      Fragment.new(primary, complement, end_for(left), end_for(right))
    end

    # The sticky end produced by a break [primary_cut, complement_cut].
    def end_for(break_pair)
      p_cut, c_cut = break_pair
      shift = signed_diff(c_cut, p_cut)
      if shift > 0
        End.new(:five_prime, strand_subseq(p_cut, p_cut + shift, false))
      elsif shift < 0
        End.new(:three_prime, strand_subseq(c_cut, c_cut - shift, false))
      else
        End.new(:blunt, '')
      end
    end

    # Sequence of the strand from cut +from+ to cut +to+.
    def strand_subseq(from, to, full_circle)
      if @circular
        length = (to - from) % seq_length
        length = seq_length if full_circle
        circular_subseq(from, length)
      else
        @sequence[from...to] || ''
      end
    end

    # +length+ bases starting at +from+, wrapping around the origin.
    def circular_subseq(from, length)
      return '' if length <= 0
      from %= seq_length
      repeats = ((from + length) / seq_length.to_f).ceil
      (@sequence * repeats)[from, length]
    end

    # Base-by-base complement of a primary-strand sequence, written
    # in the same left-to-right order (the complement strand runs
    # 3' -> 5' in this direction).
    def forward_complement(str)
      Bio::Sequence::NA.new(str).forward_complement.to_s
    end

    # Signed distance from +b+ to +a+.  For circular sequences the
    # shortest (within half the sequence length) signed distance is
    # returned.
    def signed_diff(a, b)
      diff = a - b
      if @circular && seq_length > 0
        diff %= seq_length
        diff -= seq_length if diff > seq_length / 2
      end
      diff
    end

    # Returns the length of the sequence.
    def seq_length
      @sequence.length
    end

  end # CircularAnalysis
end # RestrictionEnzyme
end # Bio
