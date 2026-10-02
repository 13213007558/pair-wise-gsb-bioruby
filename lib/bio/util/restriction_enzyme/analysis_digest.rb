#
# bio/util/restriction_enzyme/analysis_digest.rb - Explicit digestion analysis
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: Ruby License
#

require 'bio/sequence'

module Bio
class RestrictionEnzyme
class Analysis

  # Digest a linear or circular sequence using explicit caller-supplied rules.
  def self.digest(sequence, rules, options = {})
    new.digest(sequence, rules, options)
  end

  def digest(sequence, rules, options = {})
    circular = options.fetch(:circular, false)
    methylation = options[:methylation] || options[:methylations] || []
    methylations = normalize_methylations(methylation)
    rules = [rules].flatten(1).map { |rule| DigestRule.cast(rule) }
    sequence = Bio::Sequence::NA.new(sequence)

    raise ArgumentError, 'sequence is empty' if sequence.empty?
    raise ArgumentError, 'at least one enzyme rule is required' if rules.empty?

    analyzer = DigestAnalyzer.new(sequence, rules, circular, methylations)
    validate_methylation_positions(analyzer)
    analyzer.analyze
  end

  protected

  def normalize_methylations(methylations)
    Array(methylations).map do |methylation|
      case methylation
      when Methylation
        methylation
      when Integer
        Methylation.new(:+, methylation, nil)
      when Hash
        strand = methylation.fetch(:strand, :+)
        position = methylation.fetch(:position) { methylation.fetch(:pos) }
        modification = methylation[:modification] || methylation[:type]
        Methylation.new(normalize_strand(strand), position.to_i, modification)
      when Array
        if methylation.size < 2
          raise ArgumentError, "invalid methylation: #{methylation.inspect}"
        end
        Methylation.new(normalize_strand(methylation[0]),
                        methylation[1].to_i,
                        methylation[2])
      else
        raise ArgumentError, "invalid methylation: #{methylation.inspect}"
      end
    end
  end

  def normalize_strand(strand)
    case strand
    when :+, '+', :plus, 'plus' then :+
    when :-, '-', :minus, 'minus' then :-
    else raise ArgumentError, "unknown strand: #{strand.inspect}"
    end
  end

  def validate_methylation_positions(analyzer)
    analyzer.methylations.each do |methylation|
      unless methylation.position >= 0 && methylation.position < analyzer.size
        message = "methylation position #{methylation.position} is outside"
        raise ArgumentError, "#{message} the sequence"
      end
    end
  end

  # Internal coordinator for explicit rule based digestion.
  class DigestAnalyzer
    Segment = Struct.new(:strand, :start, :length, :left_cut, :right_cut)

    attr_reader :methylations, :size, :circular

    def initialize(sequence, rules, circular, methylations)
      @sequence = sequence
      @complement = sequence.forward_complement
      @rules = rules
      @circular = circular
      @methylations = methylations
      @size = sequence.length
      @result = DigestResult.new(sequence, circular)
      @cut_keys = {}
    end

    def analyze
      collect_sites
      build_fragments
      @result
    end

    protected

    def collect_sites
      @rules.each do |rule|
        scan_rule(rule, :+) if rule.scan_plus?
        scan_rule(rule, :-) if rule.scan_minus? && !rule.palindromic?
      end
      @result.cuts.sort_by! { |cut| [cut.primary_cut, cut.complement_cut] }
      @result.sticky_ends.replace(@result.cuts.map(&:sticky_end))
      @result.blocked_sites.sort_by! { |site| [site.start, site.enzyme] }
    end

    def scan_rule(rule, orientation)
      pattern = if orientation == :+
                  rule.recognition
                else
                  reverse_complement(rule.recognition)
                end
      regexp = /\A#{Bio::Sequence::NA.new(pattern).to_re.source}\z/
      positions(rule.recognition.length).each do |position|
        next unless base_slice(position, rule.recognition.length) =~ regexp

        site = site_range(position, rule.recognition.length)
        blocking = rule.blocking_methylations(@methylations,
                                             orientation,
                                             position,
                                             rule.recognition.length,
                                             @circular)
        if blocking.any?
          add_blocked_site(rule, orientation, site, blocking)
        else
          add_cuts(rule, orientation, position, site)
        end
      end
    end

    def positions(recognition_length)
      return (0...@size).to_a if @circular
      (0..@size - recognition_length).to_a
    end

    def base_slice(start, length)
      Array.new(length) { |offset| @sequence[base_index(start + offset)] }.join
    end

    def base_index(index)
      @circular ? index % @size : index
    end

    def site_range(start, length)
      if @circular
        [start, (start + length - 1) % @size]
      else
        [start, start + length - 1]
      end
    end

    def add_blocked_site(rule, orientation, site, methylations)
      @result.blocked_sites << BlockedSite.new(rule.label, orientation,
                                               site[0], site[1],
                                               methylations)
    end

    def add_cuts(rule, orientation, position, site)
      rule.physical_cuts(position, orientation).each do |primary_cut,
                                                        complement_cut|
        add_cut(rule, orientation, primary_cut, complement_cut, site)
      end
    end

    def add_cut(rule, orientation, raw_primary, raw_complement, site)
      if @circular
        primary_cut = raw_primary % @size
        complement_cut = raw_complement % @size
        displacement = circular_displacement(raw_complement - raw_primary)
      else
        validate_linear_cut(raw_primary, :primary)
        validate_linear_cut(raw_complement, :complement)
        primary_cut = raw_primary
        complement_cut = raw_complement
        displacement = raw_complement - raw_primary
      end

      key = [primary_cut, complement_cut]
      return if @cut_keys.key?(key)

      sticky_end = sticky_end(raw_primary, raw_complement, displacement)
      cut = DigestCut.new(rule.label, orientation, primary_cut,
                          complement_cut, site[0], site[1], sticky_end)
      @cut_keys[key] = cut
      @result.cuts << cut
    end

    def validate_linear_cut(cut, strand)
      return if cut >= 0 && cut <= @size - 1
      raise ArgumentError, "#{strand} cut #{cut} is outside the linear sequence"
    end

    def circular_displacement(displacement)
      displacement -= @size while displacement >= @size
      displacement += @size while displacement <= -@size
      displacement
    end

    def sticky_end(raw_primary, raw_complement, displacement)
      return StickyEnd.new(:blunt, 0, '') if displacement.zero?

      if displacement.positive?
        kind = :five_prime
        start = raw_primary + 1
      else
        kind = :three_prime
        start = raw_complement + 1
      end
      sequence = Array.new(displacement.abs) do |offset|
        @sequence[@circular ? (start + offset) % @size : start + offset]
      end.join
      StickyEnd.new(kind, displacement.abs, sequence)
    end

    def build_fragments
      if @result.cuts.empty?
        add_uncut_fragment
      else
        add_cut_fragments
      end
    end

    def add_uncut_fragment
      topology = @circular ? :circular : :linear
      @result.fragments << DigestFragment.new(@sequence.to_s,
                                             @complement.to_s,
                                             topology,
                                             0,
                                             @size - 1,
                                             0,
                                             @size - 1,
                                             nil,
                                             nil)
    end

    def add_cut_fragments
      top_cuts = @result.cuts.map { |cut| [cut.primary_cut, cut] }
      bottom_cuts = @result.cuts.map { |cut| [cut.complement_cut, cut] }
      top_segments = strand_segments(:+, top_cuts)
      bottom_segments = strand_segments(:-, bottom_cuts)
      segments = top_segments + bottom_segments
      parent = (0...segments.length).to_a
      horizontal_bases = horizontal_base_positions

      @size.times do |position|
        next if horizontal_bases.include?(position)

        top = segment_at(top_segments, position)
        bottom = segment_at(bottom_segments, position)
        union(parent, top, top_segments.length + bottom) if top && bottom
      end

      grouped_segments = group_segments(segments, parent)
      grouped_segments.each do |group|
        @result.fragments << fragment_for(group, horizontal_bases)
      end
    end

    def strand_segments(strand, boundaries)
      boundaries = boundaries.uniq { |boundary, _| boundary }.sort_by(&:first)
      return [Segment.new(strand, 0, @size, nil, nil)] if boundaries.empty?

      if @circular
        boundaries.each_with_index.map do |(boundary, cut), index|
          next_boundary, next_cut = boundaries[(index + 1) % boundaries.length]
          length = ((next_boundary - boundary - 1) % @size) + 1
          Segment.new(strand, (boundary + 1) % @size, length,
                       cut, next_cut)
        end
      else
        segments = []
        start = 0
        previous_cut = nil
        boundaries.each do |boundary, cut|
          length = boundary + 1 - start
          segments << Segment.new(strand, start, length, previous_cut, cut)
          start = boundary + 1
          previous_cut = cut
        end
        if start < @size
          segments << Segment.new(strand, start, @size - start,
                                  previous_cut, nil)
        end
        segments
      end
    end

    def horizontal_base_positions
      positions = []
      @result.cuts.each do |cut|
        displacement = cut.sticky_end.length
        next if displacement.zero?

        start = sticky_start(cut)
        displacement.times do |offset|
          positions << ((start + offset) % @size)
        end
      end
      positions
    end

    def segment_at(segments, position)
      segments.index do |segment|
        if @circular
          forward_distance(segment.start, position) < segment.length
        else
          position >= segment.start && position < segment.start + segment.length
        end
      end
    end

    def forward_distance(start, target)
      (target - start) % @size
    end

    def find(parent, index)
      parent[index] = find(parent, parent[index]) if parent[index] != index
      parent[index]
    end

    def union(parent, first, second)
      parent[find(parent, second)] = find(parent, first)
    end

    def group_segments(segments, parent)
      groups = Hash.new { |hash, key| hash[key] = [] }
      segments.each_with_index do |segment, index|
        groups[find(parent, index)] << segment
      end
      groups.values.sort_by! { |group| group_start(group) }
    end

    def fragment_for(group, horizontal_bases)
      top = group.find { |segment| segment.strand == :+ }
      bottom = group.find { |segment| segment.strand == :- }
      primary, top_start, top_end =
        strand_content(@sequence, top, horizontal_bases)
      complement, bottom_start, bottom_end =
        strand_content(@complement, bottom, horizontal_bases)
      primary, complement, top_offset, bottom_offset =
        align_strands(primary, complement, top, bottom, horizontal_bases)

      DigestFragment.new(primary,
                         complement,
                         :linear,
                         coordinate(top_start, top_offset),
                         coordinate(top_end, top_offset),
                         coordinate(bottom_start, bottom_offset),
                         coordinate(bottom_end, bottom_offset),
                         (top || bottom).left_cut,
                         (top || bottom).right_cut)
    end

    def group_start(group)
      top = group.find { |segment| segment.strand == :+ }
      bottom = group.find { |segment| segment.strand == :- }
      return top.start if top && @circular
      return bottom.start unless top
      top.start - 1
    end

    def strand_content(source, segment, _horizontal_bases)
      return ['', nil, nil] unless segment

      sequence = Array.new(segment.length) do |offset|
        index = if @circular
                  (segment.start + offset) % @size
                else
                  segment.start + offset
                end
        source[index]
      end.join
      end_position = segment.start + segment.length - 1
      end_position %= @size if @circular
      [sequence, segment.start, end_position]
    end

    def align_strands(primary, complement, top, bottom, horizontal_bases)
      return [primary, complement, 0, primary.length] if bottom.nil?
      return [primary, complement, 0, 0] if top.nil?

      paired_position = first_paired_position(top, bottom, horizontal_bases)
      top_coordinate = coordinate_distance(top, paired_position)
      bottom_coordinate = coordinate_distance(bottom, paired_position)
      relative_bottom_offset = bottom_coordinate - top_coordinate
      top_offset = [0, -relative_bottom_offset].max
      bottom_offset = relative_bottom_offset + top_offset
      width = [top.length + top_offset,
               bottom.length + bottom_offset].max
      [pad(primary, top_offset, width), pad(complement, bottom_offset, width),
       top_offset, bottom_offset]
    end

    def first_paired_position(top, bottom, horizontal_bases)
      top.length.times do |offset|
        position = if @circular
                     (top.start + offset) % @size
                   else
                     top.start + offset
                   end
        next if horizontal_bases.include?(position)
        return position if paired_position?(bottom, position)
      end
      top.start
    end

    def paired_position?(segment, position)
      if @circular
        coordinate_distance(segment, position) < segment.length
      else
        position >= segment.start &&
          position < segment.start + segment.length
      end
    end

    def sticky_start(cut)
      if cut.sticky_end.kind == :five_prime
        cut.primary_cut + 1
      else
        cut.complement_cut + 1
      end
    end

    def coordinate_distance(segment, position)
      if @circular
        forward_distance(segment.start, position)
      else
        position - segment.start
      end
    end

    def coordinate(position, _offset)
      position
    end

    def pad(string, offset, width)
      (' ' * offset) + string + (' ' * (width - string.length - offset))
    end

    def reverse_complement(sequence)
      Bio::Sequence::NA.new(sequence).reverse_complement.to_s
    end
  end

end # Analysis
end # RestrictionEnzyme
end # Bio
