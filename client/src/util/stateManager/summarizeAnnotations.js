import _ from "lodash";

/*
Build and return obs/var summary using any annotation in the schema

Summary information for each annotation, keyed by annotation name.
Value will be an object, containing summary information.

For continuous annotations (int, float, etc):
  <annotation_name>: {
    range {
      min: <number>,
      max: <number>
    }
  }

For categorical annotations (boolean, string, category):
  <annotatoin_name>: {
    options: {
      <option1>: <number>,
      ...
    },
    numOptions: <number>
  }

Summarize will be returned for BOTH obs and var annotations.

Example:
  {
    "Splice_sites_Annotated": {
      "range": {
        "min": 26,
        "max": 1075869
      }
    },
    "Selection": {
      numOptions, 6,
      "options": {
        "Astrocytes(HEPACAM)": 714,
        "Endothelial(BSC)": 123,
        "Oligodendrocytes(GC)": 294,
        "Neurons(Thy1)": 685,
        "Microglia(CD45)": 1108,
        "Unpanned": 665
      }
    }
  }

NOTE: will not summarize the required 'name' annotation, as that is
specified as unique per element.

TODO: XXX - this data structure coerces all metadata categories into a string
(ie, stores values as an Object property in the `options` field).   This looses
information (eg, type) for category types which are not strings.  Consider an
alterative data structure that does not use the object property for non-string
data types (and does not use _.countBy to summarize).
*/
function summarizeDimension(schema, annotations) {
  return _(schema)
    .filter(v => v.name !== "name")
    .keyBy("name")
    .mapValues(anno => {
      const { name, type } = anno;
      const continuous = type === "int32" || type === "float32";

      if (!continuous) {
        const categories = _.uniq(_.flatMap(annotations, name));
        const options = _.countBy(annotations, name);
        const numOptions = _.size(options);
        return {
          numOptions,
          options,
          categories
        };
      }

      if (continuous) {
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        _.forEach(annotations, obs => {
          const val = Number(obs[name]);
          min = val < min ? val : min;
          max = val > max ? val : max;
        });
        return { range: { min, max } };
      }

      throw new Error("incomprehensible schema");
    })
    .value();
}

export default function summarizeAnnotations(
  schema,
  obsAnnotations,
  varAnnotations
) {
  return {
    obs: summarizeDimension(schema.annotations.obs, obsAnnotations),
    var: summarizeDimension(schema.annotations.var, varAnnotations)
  };
}
