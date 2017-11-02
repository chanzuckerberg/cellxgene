import uri from "urijs";
import _ from "lodash";

/************************************************************************
*************************************************************************
1. This is state that will be stored in the URL query params.
1a. It is reasonably important it be kept in the same place,
because if initial load or back button etc, we'll hear about that
and manually reconstruct the state from the url, overriding whatever
we had. We could of course listen for the "url changed" type in another
reducer, if this gets messy.
1b. The data structures used here are constrained as well, as urijs
needs to be able to parse them, if arrays: https://hackernoon.com/redux-patterns-add-edit-remove-objects-in-an-array-6ee70cab2456
*************************************************************************
************************************************************************/

const selectedMetadata = (
  state = {},
  action
) => {
  switch (action.type) {
    case "url changed": {
      const {url} = action;
      const parsed = uri.parseQuery(window.location.search);
      // console.log("onURLchange sees: ", state, parsed)
      return state;
    }
    case "categorical metadata filter selected success": {
      if (!state[action.metadataField]) { /* we don't have the field, which means this is a simple insertion */
        /* {} becomes {Location: ["Tumor"]} */
        return Object.assign({}, state, {[action.metadataField]: [action.value]});
      } else { /* we do have the field, so we need to push to a copy of the array */
        const updatedField = state[action.metadataField].slice()
        updatedField.push(action.value)
        return Object.assign({}, state, {
          /* {Location: ["Tumor"]} becomes {Location: ["Tumor", "Periphery"]} */
          [action.metadataField]: updatedField
        })
      }
    }
    /*
      Seperate because it's flatter. We know in the display logic whether or not
      this the category button is selected when it's clicked (it's bold, etc),
      so we can fire a different action to save us some nesting logic here
    */
    case "categorical metadata filter deselected success": {
      /* remove the value the user just selected from the appropriate array */
      let updatedField = _.without(state[action.metadataField], action.value);
      /* if the array now looks like this {Location: []} */
      if (updatedField.length === 0) {
        return _.omit(state, action.metadataField) /* new state without Location property*/
      } else {
        return Object.assign({}, state, {
          /* {Location: ["Tumor", "Periphery"]} becomes {Location: ["Tumor"]} */
          [action.metadataField]: updatedField
        })
      }
    }
    default:
      return state;
  }

  return state;
}

export default selectedMetadata;
