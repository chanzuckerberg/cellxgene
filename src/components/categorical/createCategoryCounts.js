import * as globals from "../../globals.js"

/*
  end result:

  counts: {
    Sample.type: {
      glioblastoma: 3452
    },
    Selection: {
      foo: 234,
      bar: 21
    }
  }

*/

const categories = globals.categories;

const createCategoryCounts = (cells) => {

  /* instantiate counts obj */
  const counts = {};
  
  categories.forEach((category) => {
    counts[category] = {};
  })

  cells.forEach((cell) => {
    categories.forEach((category) => {

      /* if, for the given category (ie., categories.Location), we do not have ie., glioblastoma already, create it. Otherwise increment it. */
      if (!counts[category][cell[category]]) {
        counts[category][cell[category]] = 1;
      } else {
        counts[category][cell[category]]++
      }

    })
  });

  return counts;
}

export default createCategoryCounts;
