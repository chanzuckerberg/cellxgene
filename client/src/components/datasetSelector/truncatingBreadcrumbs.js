// Core dependencies
import { Classes, ResizeSensor } from "@blueprintjs/core";
import React, { useState } from "react";

// Characters to be used to indicate display text has been truncated
const CHAR_ELLIPSIS = "...";

// Minimum number of characters to be displayed before transitioning to a smaller state of the breadcrumbs
const MIN_VISIBLE_CHARS = 11;

// Approximate padding in pixels for each breadcrumb.
// TODO(cc) revisit - remove if we calculate actual DOM sizes rather than estimate
const ITEM_PADDING = 26;

// Approximate pixel to character ratio
const PIXELS_PER_CHAR = 6;

/*
 Individual Breadcrumb States
 ----------------------------
 F - full text
 T - truncated text
 S - indicates use of short text (eg "Collection" for collection name or "Dataset" for dataset name)
 H - hidden
 */
const STATE_FULL = "F"; // eg "Tabula Muris Senis"
const STATE_TRUNCATED = "T"; // eg "Tabula...Senis"
const STATE_SHORT_TEXT = "S"; // eg "Collection"
const STATE_HIDDEN = "H"; // --

/*
 Breadcrumbs States
 ------------------
 FFF
 FTF
 HSF
 HST
 HHS
 */
const STATES_FFF = `${STATE_FULL}${STATE_FULL}${STATE_FULL}`;
const STATES_FTF = `${STATE_FULL}${STATE_TRUNCATED}${STATE_FULL}`;
const STATES_HSF = `${STATE_HIDDEN}${STATE_SHORT_TEXT}${STATE_FULL}`;
const STATES_HST = `${STATE_HIDDEN}${STATE_SHORT_TEXT}${STATE_TRUNCATED}`;
const STATES_HHS = `${STATE_HIDDEN}${STATE_HIDDEN}${STATE_SHORT_TEXT}`;

/*
 Breadcrumb Transitions
 ----------------------
 Breadcrumbs can transition bidirectionally through states in the following order, and can also repeat individual states.
 */
const STATES = [STATES_FFF, STATES_FTF, STATES_HSF, STATES_HST, STATES_HHS];

const calculateRequiredWidth = (itemsState, items) => {
  /*
  Return the total width required to display the given items with the given states.
   */
  return items.reduce((accum, item, i) => {
    // Grab the state for this item. For example, given the state HTF, the state of the first item is H, the state of
    // the second item is T and the state of the third item is F.
    const itemState = itemsState[i];
    // Add the width (of text) corresponding to the item's state.
    if (isItemShortText(itemState)) {
      accum += item.shortTextWidth;
    } else if (isItemTruncated(itemState)) {
      accum += item.minTextWidth;
    } else if (isItemFull(itemState)) {
      accum += item.textWidth;
    }
    return accum;
  }, 0);
};

const calculateAvailableTruncatedWidth = (
  items,
  truncatedIndex,
  itemsState,
  availableWidth
) => {
  /*
  Return the width that the truncated item has available for display. That is, the available width minus the widths
  required by the other, non-truncated, items. 
   */
  // Grab the items other than the truncated item.
  const otherItems = [...items];
  otherItems.splice(truncatedIndex, 1);

  // Grab the states of the the items, other than the truncated item.
  const otherItemsState = itemsState.split("");
  otherItemsState.splice(truncatedIndex, 1);

  // Calculate the width of the other items in their corresponding states.
  const otherItemsRequiredWidth = calculateRequiredWidth(
    otherItemsState.join(""),
    otherItems
  );

  return availableWidth - otherItemsRequiredWidth;
};

const getItemsStateForAvailableWidth = (items, availableWidth) => {
  /*
  Determine the current items state (eg FFF, FTF etc) for the given available width and set of items.
   */
  for (let i = 0; i < STATES.length; i += 1) {
    const itemsState = STATES[i];
    const requiredWidth = calculateRequiredWidth(itemsState, items);
    if (availableWidth >= requiredWidth) {
      return itemsState;
    }
  }
  return STATES[STATES.length - 1]; // There's a problem, default to smallest state. TODO(cc) revisit error case here.
};

const initItems = (items) => {
  /*
  Build initial state of items, including the calculation of short text, truncated text and full text dimensions. Use
  approximation of six pixels per char. TODO(cc) revisit use of actual widths if approximation is too loose. 
   */
  return items.map((item) => {
    return {
      ...item,
      displayText: item.text, // Default display to full breadcrumb text
      minTextWidth: MIN_VISIBLE_CHARS * PIXELS_PER_CHAR + ITEM_PADDING,
      shortTextWidth: item.shortText.length * PIXELS_PER_CHAR + ITEM_PADDING,
      textWidth: item.text.length * PIXELS_PER_CHAR + ITEM_PADDING,
    };
  });
};

const isItemFull = (stateName) => {
  return stateName === STATE_FULL;
};

const isItemHidden = (stateName) => {
  return stateName === STATE_HIDDEN;
};

const isItemShortText = (stateName) => {
  return stateName === STATE_SHORT_TEXT;
};

const isItemTruncated = (stateName) => {
  return stateName === STATE_TRUNCATED;
};

const buildResizedItems = (items, availableWidth) => {
  /*
  Resize the items, either the set of visible items, or the individual item display text, to fit the given available
  width.
   */
  const itemsState = getItemsStateForAvailableWidth(items, availableWidth);
  // TODO(cc) if same state as previous and state does not contain T (eg FFF or FSF or HHS) then don't recalc here
  return updateItems(itemsState, items, availableWidth);
};

const truncate = (availableWidth, text) => {
  /*
  Return truncated text with characters removed to reduce text width to the available width. 
   */
  const visibleLength = Math.floor(availableWidth / PIXELS_PER_CHAR);
  // Determine the break indices for the "before" and "after" ellipsis text tokens
  const tokenBeforeEndIndex = Math.ceil(visibleLength / 2);
  const tokenAfterStartIndex = Math.floor(visibleLength / 2);
  // Split text at break indices and join with ellipsis
  const tokenBefore = text.substr(0, tokenBeforeEndIndex).trim();
  const tokenAfter = text.substr(text.length - tokenAfterStartIndex).trim();
  return `${tokenBefore}${CHAR_ELLIPSIS}${tokenAfter}`;
};

const updateItems = (itemsState, items, availableWidth) => {
  /*
  Update each item to match its display format to the state being transitioned to.
   */
  return items.map((item, i) => {
    const itemState = itemsState[i];
    if (isItemHidden(itemState)) {
      return {
        ...item,
        hidden: true,
      };
    }
    if (isItemShortText(itemState)) {
      return {
        ...item,
        displayText: item.shortText,
        hidden: false,
      };
    }
    if (isItemTruncated(itemState)) {
      const truncatedAvailableWidth = calculateAvailableTruncatedWidth(
        items,
        i,
        itemsState,
        availableWidth
      );
      return {
        ...item,
        displayText: truncate(truncatedAvailableWidth, item.text),
        hidden: false,
      };
    }
    return {
      ...item,
      displayText: item.text,
      hidden: false,
    };
  });
};

const TruncatingBreadcrumbs = React.memo(
  ({ breadcrumbRenderer, currentBreadcrumbRenderer, items: originalItems }) => {
    const [items, setItems] = useState(() => {
      return initItems(originalItems);
    });

    const onResize = (entries) => {
      /*
      On resize callback from ResizeSensor, save the current width of the breadcrumbs.
      */
      const availableWidth = Math.floor(entries[0].contentRect.width);
      setItems(buildResizedItems(items, availableWidth));
    };

    const renderBreadcrumb = (item, currentProp) => {
      /*
      Invoke the render callback to render the given breadcrumb. 
       */
      if (currentProp) {
        return currentBreadcrumbRenderer(item);
      }
      return breadcrumbRenderer(item);
    };

    const renderBreadcrumbs = (bcItems) => {
      /*
      Return list element/breadcrumb for each item. 
       */
      return bcItems.map((item, i) => {
        if (item.hidden) {
          return null;
        }
        // TODO(cc) possibly "clean" each item back to the format expected by BP so we can spread the Breadcrumb-specific
        // props in our render method, and so that knowledge of "displayText" vs "text"  for example, is not required by
        // parent components.
        // See datasetSelector.renderBreadcrumb for our usage, and also the following for
        // an example pattern:
        // https://github.com/palantir/blueprint/blob/826cbdf95b577c43d5fe95b99c67ee2761c853e0/packages/core/src/components/breadcrumbs/breadcrumbs.tsx#L151
        // Could possibly also have an explicit breadcrumbsProps props to neatly encapsulate and spread
        // breadcrumb-specific props, resulting in this component being a relatively transparent wrapper around
        // BP's Breadcrumbs component. For an example pattern, see `overflowListProps` on BP Breadcrumbs component.
        const currentItem = i === bcItems.length - 1;
        return <li key={item.text}>{renderBreadcrumb(item, currentItem)}</li>;
      });
    };

    return (
      <ResizeSensor onResize={onResize}>
        <ul
          className={Classes.BREADCRUMBS}
          style={{
            display: "flex",
            flexWrap: "nowrap",
            whiteSpace: "nowrap",
          }}
        >
          {renderBreadcrumbs(items)}
        </ul>
      </ResizeSensor>
    );
  }
);

export default TruncatingBreadcrumbs;
