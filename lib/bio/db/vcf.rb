# frozen_string_literal: true

require 'bio/db/vcf/header'
require 'bio/db/vcf/record'
require 'bio/db/vcf/genotype_rewriter'
require 'bio/db/vcf/field_rewriter'
require 'bio/db/vcf/allele_normalizer'
require 'bio/db/vcf/record_splitter'
require 'bio/db/vcf/normalizer'

module Bio
  # Parser and container for Variant Call Format files.
  class VCF
    META_PREFIX = '##'.freeze
    HEADER_PREFIX = '#CHROM'.freeze
    FIXED_COLUMNS = %w[#CHROM POS ID REF ALT QUAL FILTER INFO].freeze

    def initialize(str = '')
      @header = Header.new
      @records = []
      @header_seen = false
      parse(str)
    end

    attr_accessor :header, :records

    def each_record(&block)
      records.each(&block)
    end

    def normalize!(reference, sequence_start: 1, split: false)
      normalizer = Normalizer.new(reference, sequence_start: sequence_start)
      @records = records.flat_map do |record|
        normalizer.normalize(record, split: split)
      end
    end

    def normalize(reference, sequence_start: 1, split: false)
      dup.tap do |vcf|
        vcf.normalize!(reference, sequence_start: sequence_start,
                               split: split)
      end
    end

    def to_s
      lines = header.meta_lines
      lines += [header_line.join("\t") + "\n"] if @header_seen
      (lines + records.collect(&:to_s)).join
    end

    private

    def parse(str)
      str.each_line do |line|
        line = line.chomp
        next if line.empty? || (line.start_with?('#') && !meta_line?(line) &&
                               !line.start_with?(HEADER_PREFIX))

        if meta_line?(line)
          parse_meta(line)
        elsif line.start_with?(HEADER_PREFIX)
          parse_header(line)
        else
          records << Record.new(line, header)
        end
      end
    end

    def meta_line?(line)
      line.start_with?(META_PREFIX)
    end

    def parse_meta(line)
      header.meta_lines << line + "\n"
      match = line.match(/\A##(INFO|FORMAT)=<(.*)>\z/)
      return unless match

      values = parse_header_values(match[2])
      id = values['ID']
      return unless id

      target = match[1] == 'INFO' ? header.info : header.format
      target[id] = values
    end

    def parse_header_values(string)
      string.scan(/(?:\"[^\"]*\"|[^,])+/).each_with_object({}) do |part, hash|
        key, value = part.split('=', 2)
        hash[key] = value.sub(/\A\"/, '').sub(/\"\z/, '') if value
      end
    end

    def parse_header(line)
      columns = line.split("\t")
      header.sample_names = columns[9..-1].to_a
      @header_seen = true
    end

    def header_line
      return FIXED_COLUMNS if header.sample_names.empty?

      FIXED_COLUMNS + ['FORMAT'] + header.sample_names
    end
  end
end
