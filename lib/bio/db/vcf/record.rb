# frozen_string_literal: true

require 'bio/db/vcf/record_parser'
require 'bio/db/vcf/record_serializer'

module Bio
  class VCF
    # One VCF data record.
    class Record
      include RecordParser
      include RecordSerializer

      attr_accessor :chrom, :pos, :id, :ref, :alt, :qual, :filter, :info,
                    :format, :samples, :header, :normalization_status
      alias alts alt
      alias alts= alt=

      def initialize(line = nil, header = nil)
        @header = header
        @alt = []
        @info = {}
        @format = []
        @samples = {}
        @normalization_status = nil
        parse(line) if line
      end

      def initialize_copy(source)
        super
        @alt = source.alt.dup
        @info = deep_dup(source.info)
        @format = source.format.dup
        @samples = deep_dup(source.samples)
      end

      def symbolic?
        !explicit_alleles?
      end

      def explicit_alleles?
        ref.to_s.match?(/\A[A-Za-z]+\z/) &&
          alt.all? { |allele| allele.match?(/\A[A-Za-z]+\z/) }
      end

      def deep_dup(value)
        case value
        when Hash
          value.each_with_object({}) do |(key, item), hash|
            hash[key] = deep_dup(item)
          end
        when Array
          value.collect { |item| deep_dup(item) }
        when String
          value.dup
        else
          value
        end
      end
    end
  end
end
