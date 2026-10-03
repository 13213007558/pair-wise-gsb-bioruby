# frozen_string_literal: true

module Bio
  class VCF
    # Header information used while rewriting INFO and FORMAT fields.
    class Header
      attr_accessor :meta_lines, :sample_names, :info, :format

      def initialize
        @meta_lines = []
        @sample_names = []
        @info = {}
        @format = {}
      end
    end
  end
end
