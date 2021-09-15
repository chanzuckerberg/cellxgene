import React from "react";
import { Button, Tooltip, Dialog, Classes, Colors } from "@blueprintjs/core";

class AnnoDialog extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      isActive,
      text,
      title,
      instruction,
      cancelTooltipContent,
      errorMessage,
      validationError,
      annoSelect,
      annoInput,
      secondaryInstructions,
      secondaryInput,
      handleCancel,
      handleSubmit,
      primaryButtonText,
      secondaryButtonText,
      handleSecondaryButtonSubmit,
      primaryButtonProps,
    } = this.props;

    return (
      <Dialog icon="tag" title={title} isOpen={isActive} onClose={handleCancel}>
        <form
          onSubmit={(e) => {
            e.preventDefault();
          }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div style={{ marginBottom: 20 }}>
              <p>{instruction}</p>
              {annoInput || null}
              <p
                style={{
                  marginTop: 7,
                  visibility: errorMessage !== "" ? "visible" : "hidden",
                  color: Colors.ORANGE3,
                }}
              >
                {errorMessage}
              </p>
              {/* we might rename, secondary button and secondary input are not related */}
              {secondaryInstructions && (
                <p style={{ marginTop: secondaryInstructions ? 20 : 0 }}>
                  {secondaryInstructions}
                </p>
              )}
              {secondaryInput || null}
            </div>
            {annoSelect || null}
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content={cancelTooltipContent}>
                <Button onClick={handleCancel}>Cancel</Button>
              </Tooltip>
              {/* we might rename, secondary button and secondary input are not related */}
              {handleSecondaryButtonSubmit && secondaryButtonText ? (
                <Button
                  onClick={handleSecondaryButtonSubmit}
                  disabled={!text || validationError}
                  intent="none"
                  type="button"
                >
                  {secondaryButtonText}
                </Button>
              ) : null}
              <Button
                {...primaryButtonProps} // eslint-disable-line react/jsx-props-no-spreading -- Spreading props allows for modularity
                onClick={handleSubmit}
                disabled={!text || validationError}
                intent="primary"
                type="submit"
              >
                {primaryButtonText}
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    );
  }
}

export default AnnoDialog;
