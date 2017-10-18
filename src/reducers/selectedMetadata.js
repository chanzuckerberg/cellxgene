import uri from "urijs";
import _ from "lodash";

/************************************************************************
*************************************************************************
1. This is all of the state that will be stored in the URL query params.
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
    case "categorical metadata filter selected": {
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
    case "categorical metadata filter deselected": {
      return Object.assign({}, state, {

      })
    }
    default:
      return state;
  }

  return state;
}

export default selectedMetadata;
