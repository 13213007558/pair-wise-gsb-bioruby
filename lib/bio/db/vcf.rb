# frozen_string_literal: true
#
# = bio/db/vcf.rb - VCF (Variant Call Format) parser and normalizer
#
# Copyright::  Copyright (C) 2026 BioRuby project
# License::    The Ruby License
#
# == Description
#
# Parser for VCF (Variant Call Format) files, with utilities to
# normalize variants against a reference sequence (minimal
# representation and left-alignment) and to split multi-allelic
# records into bi-allelic ones.
#
# == References
#
# * VCF specification (v4.2)
#   https://samtools.github.io/hts-specs/
#

module Bio
  # Bio::VCF is a parser for VCF (Variant Call Format).
  #
  # === Examples
  #
  #   vcf = Bio::VCF.new(File.read("variants.vcf"))
  #   vcf.each do |record|
  #     record.split.each do |biallelic|
  #       puts biallelic.normalize(refseq)
  #     end
  #   end
  #
  class VCF
    include Enumerable

    # Bio::VCF::Header represents the VCF header: meta lines
    # ("##key=value"), the 8 fixed columns for the "CHROM" line, and
    # the sample names for the "FORMAT" column.
    class Header
      # Meta lines of the header as an Array of String.
      attr_reader :meta_lines

      # The column line ("#CHROM\tPOS\t...") as a String, or nil.
      attr_reader :column_line

      # Sample names as an Array of String.
      attr_reader :samples

      def initialize(meta_lines = [], column_line = nil)
        @meta_lines = meta_lines
        @column_line = column_line
        @samples = []
        parse_column_line(column_line) if column_line
        parse_meta_lines
      end

      # Returns the INFO definitions as a Hash of ID => Field.
      def infos
        @infos ||= build_fields("INFO")
      end

      # Returns the FORMAT definitions as a Hash of ID => Field.
      def formats
        @formats ||= build_fields("FORMAT")
      end

      # Returns the Number of the INFO for the given ID, or nil.
      def info_number(id)
        infos[id]&.number
      end

      # Returns the Number of the FORMAT for the given ID, or nil.
      def format_number(id)
        formats[id]&.number
      end

      # Returns the parsed header as a String.
      def to_s
        lines = @meta_lines
        lines = [@column_line] if @column_line.nil?
        lines.join("\n")
      end

      private

      def parse_column_line(line)
        cols = line.chomp.split("\t")
        @samples = cols[9..] || []
      end

      def parse_meta_lines
        @infos = nil
        @formats = nil
      end

      def build_fields(kind)
        fields = {}
        @meta_lines.each do |line|
          next unless line.start_with?("##")
          id = parse_meta_line(line)
          next unless id
          fields[id] = Field.new(id, nil, nil, nil)
        end
        fields
      end

      def parse_meta_line(line)
        return nil unless line =~ /\A##([^=]+)=<(.+)>\z/m
        id = Regexp.last_match(1)
        fields = {}
        inner = Regexp.last_match(2)
        inner.split(/,(?=(?:[^"]|"[^"]*")*$)/).each do |pair|
          k, v = pair.split("=", 2)
          fields[k] = v
        end
        [id, fields]
      end
    end

    # Bio::VCF::Record represents a single VCF record (one line in
    # the two-line two-column tabular
    # for each alternate allele.
    class Record
      # Creates a new Record object.
      def initialize(chrom, pos, id, ref, alt, qual, filter, info, format, samples)
        @chrom = chrom
        @pos = pos
        @id = id
        @ref = ref
        @alt = alt
        @qual = qual
        @filter = filter
        @info = info
        @format = format
        @samples = samples
      end

      # Chromosome name.
      attr_accessor :chrom

      # 1-based position of the variant on the chromosome.
      attr_accessor :pos

      # Identifier of the variant.
      attr_accessor :id

      # Reference allele.
      attr_accessor :ref

      # Alternate alleles (Array of String).
      attr_accessor :alt

      # Quality score.
      attr_accessor :qual

      # Filter status.
      attr_accessor :filter

      # INFO field as a Hash of String => String or Array of String.
      attr_accessor :info

      # FORMAT field as a String.
      attr_accessor :format

      # Sample data as an Array of Hash for each sample.
      attr_accessor :samples

      # Returns true if the record is a multi-allelic record.
      def multi_allelic?
        alt.size > 1
      end

      # Returns true if the record is a bi-allelic record.
      def bi_allelic?
        alt.size == 1
      end

      # Returns true if the record contains a symbolic allele
      # (e.g. "<DEL>", "<INS>", "<DUP>").
      def symbolic?
        alt.any? { |a| a.include?("<") || a.include?("[") || a.include?("]") }
      end

      # Returns true if the record contains alt alleles are all explicit
      # sequence alleles (A, C, G, T, N, or the empty string for
      # deletions).
      def explicit?
        alt.all? { |a| a =~ /\A[ACGTNacgtn]+\z/ || a == "" }
      end

      # Returns true if the record contains only explicit sequence
      alleles and no symbolic alleles.
      def processable?
        !symbolic? && explicit?
      end

      # Splits a multi-allelic record into bi-allelic records.
      # Returns an Array of Record objects.
      # If the record is already bi-allelic, returns [self].
      # If the record contains symbolic alleles, returns [self]
      # unchanged (symbolic alleles are not processed).
      def split
        return [self] unless multi_allelic?
        return [self] unless processable?

        alt.flat_map.with_index do |a, i|
          rec = dup_record
          rec.alt = [a]
          rec.info = split_info(info, i)
          rec.samples = split_samples(samples, i)
          rec
        end
      end

      # Normalizes the record against the reference sequence.
      # Returns a new Record object with minimal representation
      # (left-aligned and trimmed).
      # If the record contains symbolic alleles, returns self
      # unchanged (symbolic alleles are not processed).
      def normalize(refseq)
        return self unless processable?

        ref = refseq.is_a?(String) ? refseq[pos - 1, 1] : nil
        return self unless ref

        normalized = dup_record
        normalized.normalize!(refseq)
        normalized
      end

      # Normalizes the record in place.
      # Returns self.
      def normalize!(refseq)
        return self unless processable?

        # Trim common suffix
        trim_suffix!(refseq)
        # Trim common prefix
        trim_prefix!(refseq)
        # Left-align
        left_align!(refseq)
        self
      end

      # Returns the normalized record as a String in VCF format.
      def to_s
        if format
          [chrom, pos, id, ref, alt.join(","), qual, filter, info_str, format, *samples_str].join("\t")
        else
          [chrom, pos, id, ref, alt.join(","), qual, filter, info_str].join("\t")
        end
      end

      # Returns the INFO field as a String.
      def info_str
        if info.is_a?(Hash)
          info.map { |k, v| v.is_a?(Array) ? "#{k}=#{v.join(",")}" : "#{k}=#{v}" }.join(";")
        else
          info
        end
      end

      # Returns the samples as an Array of String.
      def samples_str
        samples.map do |s|
          if s.is_a?(Hash)
            s.map { |k, v| v.is_a?(Array) ? "#{k}=#{v.join(",")}" : "#{k}=#{v}" }.join(":")
          else
            s
          end
        end
      end

      private

      def dup_record
        rec = dup
        rec.instance_variable_set(:@alt, alt.dup) if alt.is_a?(Array)
        rec.instance_variable_set(:@info, info.dup) if info.is_a?(Hash)
        rec.instance_variable_set(:@samples, samples.dup) if samples.is_a?(Array)
        rec
      end

      def split_info(info, index)
        return info unless info.is_a?(Hash)

        info.each_with_object({}) do |(k, v), h|
          h[k] = if v.is_a?(Array)
                   v[index] || v.first
                 else
                   v
                 end
        end
      end

      def split_samples(samples, index)
        return samples unless samples.is_a?(Array)

        samples.map do |s|
          if s.is_a?(Hash)
            s.each_with_object({}) do |(k, v), h|
              h[k] = if v.is_a?(Array)
                       v[index] || v.first
                     else
                       v
                     end
            end
          else
            s
          end
        end
      end

      def trim_suffix!(refseq)
        return unless refseq.is_a?(String)

        ref = refseq[pos - 1, 1]
        return unless ref

        # Trim common suffix from all alleles
        while ref.length > 1 && alt.all? { |a| a.length > 1 } && alt.all? { |a| a[-1] == ref[-1] }
          ref = ref[0..-2]
          alt.each { |a| a.chop! }
        end
      end

      def trim_prefix!(refseq)
        return unless refseq.is_a?(String)

        ref = refseq[pos - 1, 1]
        return unless ref

        # Trim common prefix from all alleles
        while ref.length > 1 && alt.all? { |a| a.length > 1 } && alt.all? { |a| a[0] == ref[0] }
          ref = ref[1..]
          alt.each { |a| a[0] = "" }
        end
      end

      def left_align!(refseq)
        return unless refseq.is_a?(String)

        ref = refseq[pos - 1, 1]
        return unless ref

        # Left-align: if the first base of REF is the same as the
        # first base of ALT, and the last base of REF is the same as
        # the last base of ALT, we can shift left.
        while ref.length > 1 && alt.all? { |a| a.length > 1 } &&
              alt.all? { |a| a[0] == ref[0] } && alt.all? { |a| a[-1] == ref[-1] }
          ref = ref[1..]
          alt.each { |a| a[0] = "" }
        end
      end
    end

    # Creates a new VCF object.
    # If +str+ is given, parses it.
    def initialize(str = nil)
      @header = nil
      @records = []
      if str
        parse(str)
      end
    end

    # Parses the given string and returns self.
    def parse(str)
      @header = nil
      @records = []
      lines = str.split(/\r?\n/)
      header_lines = []
      record_lines = []
      lines.each do |line|
        if line.start_with?("#")
          header_lines << line
        else
          record_lines << line
        end
      end
      @header = Header.new
      @header.parse(str) if header_lines.any?
      @records = record_lines.map { |l| parse_record(l) }
      @records
    end

    # Returns the header as a Header object.
    attr_reader :header

    # Returns the records as an Array of Record objects.
    attr_reader :records

    # Iterates over the records.
    def each(&block)
      @records.each(&block)
    end

    # Returns the VCF as a String.
    def to_s
      if header
        header.to_s + "\n" + records.map(&:to_s).join("\n")
      else
        records.map(&:to_s).join("\n")
      end
    end

    private

    def parse_record(line)
      fields = line.split("\t")
      chrom = fields[0]
      pos = fields[1].to_i
      id = fields[2]
      ref = fields[3]
      alt = fields[4].split(",")
      qual = fields[5]
      filter = fields[6]
      info = fields[7]
      format = fields[8]
      samples = fields[9..] || []
      Record.new(chrom, pos, id, ref, alt, qual, filter, info, format, samples)
    end
  end
end
